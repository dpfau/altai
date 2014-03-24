function ROI = detectROIs(params)

numROI = int32(0);
datadim = length(params.sz);

patchrng = arrayfun(@(x)1:x,params.roiSz,'UniformOutput',0);

ROIShapes = zeros([params.roiSz,params.maxROI]); % Initialize the whole array. Expanding as we go is slow and dumb.
ROIPrecs  = zeros([params.roiSz,params.maxROI]); % the precision of each pixel in the ROI.
ROIOffset = zeros(datadim,params.maxROI,'int32'); % Location of the ROIs. Fixed from the start, integer precision
ROICenter = zeros(datadim,params.maxROI); % Slightly different from ROIOffset. Double precision, updated online, used to decide if two ROI are close enough to merge.
ROIPower  = zeros(1,params.maxROI); % sum of squared firing rates over all ROIs
ROITimes = sparse(1000,1e5); % index all the times at which a given ROI appears
OutOfBounds = 0; % track the number of ROIs we toss out (should be negligible)

vec = @(x)x(:);
ROIRng = @(x) arrayfun(@(x,y)x-int32(floor(y/2))+(0:int32(y)-1),x,params.roiSz','UniformOutput',0);

for t = params.tRng
    tic
	watersheds = zeros(params.sz,'int32');
    data = params.loadframe(t);
    fprintf('%d: ',t);
    if gpuDeviceCount
        gpuDataBlur = blur(gpuArray(data), params.sig);
        fprintf('B');

        gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
        regmax = gather(gpuRegmax);
    else
        dataBlur = blur(data, params.sig);
        fprintf('B');

        regmax = int32(find(myregionalmax(dataBlur-params.thresh)));
    end
    numRegmax = length(regmax);
    
    if datadim == 2
        [xRegmax, yRegmax] = ind2sub(params.sz,regmax);
        regmaxSub = [xRegmax,yRegmax];
    else
        [xRegmax, yRegmax, zRegmax] = ind2sub(params.sz,regmax);
        regmaxSub = [xRegmax,yRegmax,zRegmax];
    end
    
    
    inBounds = ~any( regmaxSub - int32(floor(params.roiSz(ones(numRegmax,1),:)/2))<=0 | ...
                     regmaxSub + int32(floor(params.roiSz(ones(numRegmax,1),:)/2))>params.sz(ones(numRegmax,1),:), 2 );

    regmax = regmax(inBounds);
    % Update gpuRegmax?
    gpuRegmax = regmax;
    xRegmax = xRegmax(inBounds); yRegmax = yRegmax(inBounds); 
    
    if datadim == 3
        zRegmax = zRegmax(inBounds);
    end
    regmaxSub = double(regmaxSub(inBounds,:));

    if gpuDeviceCount        
        intensity = gather(gpuDataBlur(gpuRegmax));
    else
        intensity = dataBlur(regmax);
    end
    OutOfBounds = OutOfBounds + sum(~inBounds);
    fprintf('M');

    if params.watershedFlag
        if gpuDeviceCount
            fastwatershed(gather(gpuDataBlur-params.thresh), watersheds, regmax);
            clear gpuDataBlur
            watersheds = gpuArray(watersheds);
        else
            fastwatershed(dataBlur-params.thresh, watersheds, regmax);
        end
    else
        if gpuDeviceCount
            dataBlur = gather(gpuDataBlur);
        end
        watersheds = double(watershed(-dataBlur)) .* (dataBlur > params.thresh);
    end
    fprintf('W');
    
    % DEBUG: Set the threshold based on watershed.
    h = figure();
    subplot(2,1,1);
    imagesc(gather(gpuDataBlur));
    colorbar();
    colormap gray;
    axis image;
    subplot(2,1,2);
    imagesc(watersheds);
    hold on
    scatter(yRegmax, xRegmax, 'ro');
    colorbar();
    axis image;
    saveas(h, sprintf('watershed_%d.png', t));
    close(h);

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
        [nearestNeighbors, nnDistance] = knnsearch((ROICenter(:,1:numROI))'/diag(params.sig), regmaxSub/diag(params.sig), 'K', 1);
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
                region(rng{:}) = region(rng{:}) & ROIShapes(patchrng{:},idx(i));
                if nnz(region)>30
                    vars(i) = tryGather( var(residual(region)) );
                end
            end
            foo = pinv([new_rates+new_rates.^2./old_rates, 1+new_rates.^2./old_rates.^2])*vars;
            params.var = foo(1);
            params.varSlope = foo(2);
            for i = 1:numROI
                rng = ROIRng(ROIOffset(:,i));
                ROIPrecs(patchrng{:},i) = old_intensity(i)^2 .* logical(ROIShapes(patchrng{:},i)) / (params.var + params.varSlope*old_intensity(i));
            end
        end

        % Scale residual to be unit norm
        residVar = params.var*ones(size(residual));
        for i = 1:numROI
            rng = ROIRng(ROIOffset(:,i));
            residVarUpdate = params.var * rates(i) * logical(ROIPrecs(patchrng{:},i)) + rates(i)^2./ROIPrecs(patchrng{:},i);
            residVarUpdate(~logical(ROIPrecs(patchrng{:},i))) = 0;
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
        allNeighbors = rangesearch((ROICenter(:,1:numROI))'/diag(params.roiSz), regmaxSub/diag(params.roiSz), 1, 'NSMethod', 'kdtree', 'Distance', 'chebychev');
        numChi2 = 0; % number of regional maxima merged into existing ROIs because they passed a Chi^2 test
        for i = 1:length(regmax)
            if assignment(i) == 0
                % Get region within watershed that overlaps other ROIs
                if gpuDeviceCount
                    region = gpuArray.false(size(watersheds));
                else
                    region = false(size(watersheds));
                end
                for j = 1:length(allNeighbors{i})
                    idx = allNeighbors{i}(j);
                    rng = ROIRng(ROIOffset(:,idx));
                    region(rng{:}) = region(rng{:}) | ROIShapes(patchrng{:},idx);
                end
                region = region & watersheds == i;

                % See if residual passes Chi^2 test
                error1 = tryGather(norm(residual(region))^2); % scale of the residual in the overlap between ROI and watershed
                distance = ROICenter(:,nearestNeighbors(i))-double(regmaxSub(i,:)');
                error2 = distance'/diag(params.sig)/params.minDist^2*distance; % scaled distance between the nearest neighbor ROI and the regional maximum
                if (1 - chi2cdf(error1,tryGather(nnz(region)))) * (1 - chi2cdf(error2,3)) > params.pval
                    % Assign regional maximum to ROI with greatest power
                    pow = zeros(length(allNeighbors{i}),1);
                    for j = 1:length(allNeighbors{i})
                        idx = allNeighbors{i}(j);
                        rng = ROIRng(ROIOffset(:,idx));
                        ROIShape = ROIShapes(patchrng{:},idx);
                        pow(j) = norm(rates(idx)*vec(ROIShape(region(rng{:}))));
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
                if datadim == 2
                    [ROICenter(1,numROI), ROICenter(2,numROI)] = ind2sub(params.sz,regmax(i));
                elseif datadim == 3
                    [ROICenter(1,numROI), ROICenter(2,numROI), ROICenter(3,numROI)] = ind2sub(params.sz,regmax(i));
                end
                ROIOffset(:,numROI) = int32(ROICenter(:,numROI));
                rng = ROIRng(ROIOffset(:,numROI));
                ROIShapes(patchrng{:},numROI) = tryGather( residual(rng{:}) .* sqrt(residVar(rng{:})) .* (watersheds(rng{:})==i) / intensity(i) ); % remember to re-scale the residual back to what it originally was
                ROIPrecs(patchrng{:},numROI) = tryGather( intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i)) );
                ROIPower(numROI) = intensity(i)^2;
                ROITimes(t,numROI) = 1;
            else % merge ROI
                j = assignment(i);
                rng = ROIRng(ROIOffset(:,j));
                ROIShapes(patchrng{:},j) = tryGather( (ROIPrecs(patchrng{:},j) .* ROIShapes(patchrng{:},j) + ...
                    intensity(i)*(watersheds(rng{:})==i)/(params.var + params.varSlope*intensity(i)) .* data(rng{:})) ./ ... % this line right here might be why sometimes we get multiple ROIs mixed together. And why aren't we updating all ROIs?
                (ROIPrecs(patchrng{:},j) + intensity(i)^2*(watersheds(rng{:})==i)/(params.var + params.varSlope*intensity(i))) );
                ROICenter(:,j) = (ROIPower(j) * ROICenter(:,j) + intensity(i).^2 * double(regmaxSub(i,:)'))...
                    /(ROIPower(j) + intensity(i).^2);
                ROIPower(j) = ROIPower(j) + intensity(i)^2;
                ROIPrecs(patchrng{:},j) = tryGather( ROIPrecs(patchrng{:},j) + intensity(i)^2*(watersheds(rng{:})==i)/(params.var + params.varSlope*intensity(i)) );
                ROITimes(t,j) = 1;
            end
        end
        ROIShapes(isnan(ROIShapes)) = 0;
    else
        numROI = numROI + length(regmax);
        ROITimes(t,1:numROI) = 1;
        for i = 1:numROI
            if datadim == 2
                [ROICenter(1,i), ROICenter(2,i)] = ind2sub(params.sz,regmax(i));
            else
                [ROICenter(1,i), ROICenter(2,i), ROICenter(3,i)] = ind2sub(params.sz,regmax(i));
            end
            ROIOffset(:,i) = int32(ROICenter(:,i));
            rng = ROIRng(ROIOffset(:,i));
            ROIShapes(patchrng{:},i) = tryGather( data(rng{:}) .* (watersheds(rng{:})==i) / intensity(i) );
            if ~params.autoVar
                ROIPrecs(patchrng{:},i) = tryGather( intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i)) );
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
ROIShapes = ROIShapes(patchrng{:},1:numROI);
ROIPrecs  = ROIPrecs(patchrng{:},1:numROI);
ROIOffset = ROIOffset(:,1:numROI);
ROICenter = ROICenter(:,1:numROI);
ROIPower  = ROIPower(1:numROI);
ROI = makeROIStruct(numROI, ROIShapes, ROIPrecs, ROIOffset, ROICenter, ROIPower, ROITimes);

function ROI = makeROIStruct(numROI,ROIShapes,ROIPrecs,ROIOffset,ROICenter,ROIPower,ROITimes)

% Create a struct array where each ROI is one struct, instead of keeping all the data in separate
% fields. This is really an intermediate step, that should precede the whole function being rewritten
% to work with structs from the start

roiSz = size(ROIShapes);
roiSz = roiSz(1:end-1);
patchrng = arrayfun(@(x)1:x,roiSz,'UniformOutput',0);

ROI = [];
for i = 1:numROI
    ROI(i).shape = ROIShapes(patchrng{:},i);
    ROI(i).prec  = ROIPrecs(patchrng{:},i);
    ROI(i).offset = ROIOffset(:,i);
    ROI(i).center = ROICenter(:,i);
    ROI(i).power = ROIPower(i);
    ROI(i).times = ROITimes(:,i);
end
