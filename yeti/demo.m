clear params
params.sig = 5;
params.dz = 12.5;
params.thresh = 0.01;
params.roiSz = [40,40,3]; % Max size of an ROI
params.sz = [1024,2048,43];
params.minDist = 3; % minimum distance below which two ROI are considered the same
params.autoVar = true; % automatically compute baseline variance and linear dependence on firing rate
params.maxROI = 1e5;
params.pval = 1e-14; % very strict.
params.posVar = 5; % Variance in the ROI position. Should set higher than the actual variance, to be robust to noise, but still low enough that it stops outrageous things.

numROI = int32(0);
ROIShapes = zeros([params.roiSz,params.maxROI]); % Initialize the whole sparse array. Expanding as we go is slow and dumb.
ROIPrecs  = zeros([params.roiSz,params.maxROI]); % the precision of each pixel in the ROI.
ROIOffset = zeros(3,params.maxROI,'int32'); % Location of the ROIs. Fixed from the start, integer precision
ROICenter = zeros(3,params.maxROI); % Slightly different from ROIOffset. Double precision, updated online, used to decide if two ROI are close enough to merge.
ROIPower  = zeros(1,params.maxROI); % sum of squared firing rates over all ROIs
OutOfBounds = 0; % track the number of ROIs we toss out (should be negligible)

vec = @(x)x(:);
ROIRng = @(x) arrayfun(@(x,y)x-int32(floor(y/2))+(0:int32(y)-1),x,params.roiSz','UniformOutput',0);

if ~exist('tRng','var'), tRng = [25:1000,1:24]; end
for t = tRng
    tic
	watersheds = zeros(params.sz,'int32');
    if gpuDeviceCount
        data = padarray(loadframe(t),[0,0,1]); % pad the third dimenions, in case some ROIs bump against the edge (more likely than in the other two dimensions)
        fprintf('%d: ',t);
        gpuDataBlur = blur(gpuArray(data), [params.sig, params.sig, params.sig/params.dz]);
        fprintf('B');

        gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
        regmax = gather(gpuRegmax);
    else
        data = padarray(loadframe(t,'/Users/pfau/Documents/Research/ROI/Light Sheet Data/binary_data'),[0,0,1]);
        fprintf('%d: ',t);
        dataBlur = blur(data, [params.sig, params.sig, params.sig/params.dz]);
        fprintf('B');

        regmax = int32(find(myregionalmax(dataBlur-params.thresh)));
    end
    numRegmax = length(regmax);
    [xRegmax, yRegmax, zRegmax] = ind2sub(params.sz,regmax);
    regmaxSub = [xRegmax,yRegmax,zRegmax];
    inBounds = ~any( regmaxSub + int32(floor(params.roiSz(ones(numRegmax,1),:)/2))<0 | regmaxSub + int32(floor(params.roiSz(ones(numRegmax,1),:)/2))>params.sz(ones(numRegmax,1),:), 2 );
    regmax = regmax(inBounds);
    xRegmax = xRegmax(inBounds); yRegmax = yRegmax(inBounds); zRegmax = zRegmax(inBounds);
    regmaxSub = regmaxSub(inBounds,:);
    if gpuDeviceCount
        intensity = gather(gpuDataBlur(gpuRegmax));
    else
        intensity = dataBlur(regmax);
    end
    OutOfBounds = OutOfBounds + sum(~inBounds);
    fprintf('M');

    if gpuDeviceCount
        fastwatershed(gather(gpuDataBlur-params.thresh), watersheds, regmax);
        clear gpuDataBlur
        watersheds = gpuArray(watersheds);
    else
        fastwatershed(dataBlur-params.thresh, watersheds, regmax);
    end
    fprintf('W');

    % should probably drop this section into its own function
    if numROI > 0
        assignment = zeros(length(regmax),1); % index of ROI to assign regional maximum to, or 0 if it's a new ROI
        % Compute residual
        [rates,residual] = ratesFromFrame(data,ROIShapes,ROIOffset,numROI);
        if gpuDeviceCount
            residual = gpuArray(residual);
        end

        % Compute nearest neighbors, if regional maxima are close enough to ROI centers, merge them together
        warning('off','stats:KDTreeSearcher:knnsearch:DataConversion');
        [nearestNeighbors, nnDistance] = knnsearch((diag([1,1,params.dz])*ROICenter(:,1:numROI))', [xRegmax,yRegmax,zRegmax*params.dz], 'K', 1);
        assignment(nnDistance < params.minDist) = nearestNeighbors(nnDistance < params.minDist);
        numNeighbors = nnz(nnDistance < params.minDist);
        if params.autoVar && ~isfield(params,'varSlope') % estimate the dependence of the variance on the firing rate
            idx = find(nnDistance < params.minDist);
            vars = zeros(numNeighbors,1);
            old_rates = old_intensity(assignment(nnDistance < params.minDist));
            new_rates = intensity(nnDistance < params.minDist);
            for i = 1:numNeighbors
                region = watersheds==idx(i);
                rng = ROIRng(ROIOffset(:,assignment(idx(i))));
                region(rng{:}) = region(rng{:}) & ROIShapes(:,:,:,idx(i));
                vars(i) = tryGather( var(residual(region)) );
            end
            foo = pinv([new_rates+new_rates.^2./old_rates, 1+new_rates.^2./old_rates.^2])*vars;
            params.var = foo(1);
            params.varSlope = foo(2);
            for i = 1:numROI
                rng = ROIRng(ROIOffset(:,i));
                ROIPrecs(:,:,:,i) = tryGather( old_intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*old_intensity(i)) );
            end
        end

        % Scale residual to be unit norm
        residVar = params.var*ones(size(residual));
        for i = 1:numROI
            rng = ROIRng(ROIOffset(:,i));
            residVarUpdate = params.var * rates(i) * logical(ROIPrecs(:,:,:,i)) + rates(i)^2./ROIPrecs(:,:,:,i);
            residVarUpdate(~logical(ROIPrecs(:,:,:,i))) = 0;
            try
                residVar(rng{:}) = residVar(rng{:}) + residVarUpdate;
            catch e
                disp(e);
            end
        end
        residual = residual ./ sqrt(residVar);
        fprintf('R');
        
        % Get index of ROIs that overlap regional maxima
        warning('off','stats:KDTreeSearcher:rangesearch:DataConversion'); % don't need to hear about my conversions
        allNeighbors = rangesearch((diag([1,1,params.dz])*ROICenter(:,1:numROI))',[xRegmax,yRegmax,zRegmax*params.dz], max(params.roiSz .* [1,1,params.dz]), 'NSMethod', 'kdtree', 'Distance', 'chebychev');
        numChi2 = 0; % number of regional maxima merged into existing ROIs because they passed a Chi^2 test
        for i = 1:length(regmax)
            if assignment(i) == 0
                % Get region within watershed that overlaps other ROIs
                rng = ROIRng([xRegmax(i);yRegmax(i);zRegmax(i)]);
                if gpuDeviceCount
                    region = gpuArray.false(size(watersheds));
                else
                    region = false(size(watersheds));
                end
                for j = 1:length(allNeighbors{i})
                    region(rng{:}) = region(rng{:}) | ROIShapes(:,:,:,allNeighbors{i}(j));
                end
                region = region & watersheds == i;

                % See if residual passes Chi^2 test
                error1 = tryGather(norm(residual(region))^2); % scale of the residual in the overlap between ROI and watershed
                distance = ROICenter(:,nearestNeighbors(i))-double([xRegmax(i);yRegmax(i);zRegmax(i)]);
                error2 = distance'/diag([params.minDist^2, params.minDist^2, (params.minDist/params.dz)^2])*distance; % scaled distance between the nearest neighbor ROI and the regional maximum
                if (1 - chi2cdf(error1,tryGather(nnz(region)))) * (1 - chi2cdf(error2,3)) > params.pval
                    % Assign regional maximum to ROI with greatest power
                    pow = zeros(length(allNeighbors{i}),1);
                    for ii = 1:length(allNeighbors{i})
                        idx = allNeighbors{i}(ii);
                        rng = ROIRng(ROIOffset(:,idx));
                        ROIShape = ROIShapes(:,:,:,idx);
                        pow(ii) = norm(rates(idx)*vec(ROIShape(region(rng{:}))));
                    end
                    [~,jj] = max(pow);
                    assignment(i) = allNeighbors{i}(jj);
                    numChi2 = numChi2 + 1;
                end
            end
        end

        for i = 1:length(regmax)
            if assignment(i) == 0 % create new ROI
                numROI = numROI + 1;
                if numROI > params.maxROI, pause; end

                [ROICenter(1,numROI), ROICenter(2,numROI), ROICenter(3,numROI)] = ind2sub(params.sz,regmax(i));
                ROIOffset(:,numROI) = int32(ROICenter(:,numROI));
                rng = ROIRng(ROIOffset(:,numROI));
                ROIShapes(:,:,:,numROI) = tryGather( residual(rng{:}) .* (watersheds(rng{:})==i) / intensity(i) );
                ROIPrecs(:,:,:,numROI) = tryGather( intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i)) );
                ROIPower(numROI) = intensity(i)^2;
            else % merge ROI
                j = assignment(i);
                rng = ROIRng(ROIOffset(:,j));
                ROIShapes(:,:,:,j) = tryGather( (ROIPrecs(:,:,:,j) .* ROIShapes(:,:,:,j) + ...
                    intensity(i)*(watersheds(rng{:})==i)/(params.var + params.varSlope*intensity(i)) .* data(rng{:})) ./ ... % this line right here might be why sometimes we get multiple ROIs mixed together. And why aren't we updating all ROIs?
                (ROIPrecs(:,:,:,j) + intensity(i)^2*(watersheds(rng{:})==i)/(params.var + params.varSlope*intensity(i))) );
                ROICenter(:,j) = (ROIPower(j) * ROICenter(:,j) + intensity(i).^2 * double([xRegmax(i); yRegmax(i); zRegmax(i)]))...
                    /(ROIPower(j) + intensity(i).^2);
                ROIPower(j) = ROIPower(j) + intensity(i)^2;
                ROIPrecs(:,:,:,j) = tryGather( ROIPrecs(:,:,:,j) + intensity(i)^2*(watersheds(rng{:})==i)/(params.var + params.varSlope*intensity(i)) );
            end
        end
        ROIShapes(isnan(ROIShapes)) = 0;
    else
        numROI = numROI + length(regmax);
        for i = 1:numROI
            [ROICenter(1,i), ROICenter(2,i), ROICenter(3,i)] = ind2sub(params.sz,regmax(i));
            ROIOffset(:,i) = int32(ROICenter(:,i));
            rng = ROIRng(ROIOffset(:,i));
            ROIShapes(:,:,:,i) = tryGather( data(rng{:}) .* (watersheds(rng{:})==i) / intensity(i) );
            if ~params.autoVar
                ROIPrecs(:,:,:,i) = tryGather( intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i)) );
            end
        end
        ROIPower(1:numROI) = intensity.^2;
        numNeighbors = 0;
        numChi2 = 0;
        % Automatically set the baseline variance
        if params.autoVar
            old_intensity = intensity; % save the intensity from the first frame to compute the dependence of variance on firing rate
        end
    end
    fprintf('\t Found %d regional maxima: %d neighbors, %d pass Chi^2, %d new, %d total ROI, %gs\n', length(regmax), numNeighbors, numChi2, length(regmax)-numNeighbors-numChi2, numROI, toc);
end
ROIShapes = ROIShapes(:,:,:,1:numROI);
ROIPrecs  = ROIPrecs(:,:,:,1:numROI);
ROIOffset = ROIOffset(:,1:numROI);
ROICenter = ROICenter(:,1:numROI);
ROIPower  = ROIPower(1:numROI);
save ROI_results numROI ROIShapes ROIPrecs ROIOffset ROICenter ROIPower OutOfBounds
