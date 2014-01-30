params.sig = 5;
params.dz = 12.5;
params.thresh = 0.011;
params.roiSz = [40,40,3]; % Max size of an ROI
params.sz = [1472,2048,41] + 2*floor(params.roiSz/2); % don't forget the padding!
params.minDist = 3; % minimum distance below which two ROI are considered the same
params.var = 1; % baseline variance
params.varSlope = 0.001; % slope of the variance as a function of the intensity
params.maxROI = 1e5;
params.pval = 1e-14; % very strict.

numROI = int32(0);
ROIShapes = zeros([params.roiSz,params.maxROI]); % Initialize the whole sparse array. Expanding as we go is slow and dumb.
ROIPrecs  = zeros([params.roiSz,params.maxROI]); % the precision of each pixel in the ROI.
ROIOffset = zeros(3,params.maxROI,'int32'); % Location of the ROIs. Fixed from the start, integer precision
ROICenter = zeros(3,params.maxROI); % Slightly different from ROIOffset. Double precision, updated online, used to decide if two ROI are close enough to merge.
ROIPower  = zeros(1,params.maxROI); % sum of squared firing rates over all ROIs

vec = @(x)x(:);
ROIRng = @(x) arrayfun(@(x,y)x-floor(int32(y)/2)+(0:int32(y)-1),x,params.roiSz','UniformOutput',0);

for t = 100:110
    tic
	watersheds = zeros(params.sz,'int32');
    data = padarray(loadframe(t),floor(params.roiSz/2)); % Prevents ROIs from going over the edge of an image.
    fprintf('%d: ',t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('B');
    clear gpuData; % save space on the GPU
    gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
    regmax = gather(gpuRegmax);
    intensity = gather(gpuDataBlur(gpuRegmax));
    fprintf('M');
    fastwatershed(gather(gpuDataBlur-params.thresh), watersheds, regmax);
    fprintf('W');
    % should probably drop this section into its own function
    if numROI > 0
        assignment = zeros(length(regmax),1); % index of ROI to assign regional maximum to, or 0 if it's a new ROI
        % Compute residual
        residual = double(data); fprintf('.');
        rates = ratesfromframe(residual, ROIShapes, ROIOffset, numROI); fprintf('.');
        getResidual(data,residual,ROIShapes,ROIOffset,rates); fprintf('.');
        % Scale residual to be unit norm
        residVar = params.var*ones(size(residual));
        for i = 1:numROI
            rng = ROIRng(ROIOffset(:,i));
            residVarUpdate = params.var * rates(i) * logical(ROIPrec(:,:,:,i)) + rates(i)^2./ROIPrec(:,:,:,j);
            residVarUpdate(~logical(ROIPrec(:,:,:,i))) = 0;
            residVar(rng{:}) = residVar(rng{:}) + residVarUpdate;
        end
        residual = residual ./ sqrt(residVar);
        fprintf('R');

        % Compute nearest neighbors, if regional maxima are close enough to ROI centers, merge them together
        [xRegmax, yRegmax, zRegmax] = ind2sub(params.sz,regmax);
        warning('off','stats:KDTreeSearcher:knnsearch:DataConversion');
        [nearestNeighbors, nnDistance] = knnsearch((diag([1,1,params.dz])*ROICenter(:,1:numROI))', [xRegmax,yRegmax,zRegmax*params.dz], 'K', 1);
        assignment(nnDistance < params.minDist) = nearestNeighbors(nnDistance < params.minDist);
        numNeighbors = nnz(nnDistance < params.minDist);

        % Get index of ROIs that overlap regional maxima
        warning('off','stats:KDTreeSearcher:rangesearch:DataConversion'); % don't need to hear about my conversions
        allNeighbors = rangesearch((diag([1,1,params.dz])*ROICenter(:,1:numROI))',[xRegmax,yRegmax,zRegmax*params.dz], max(params.sz .* [1,1,params.dz]), 'NSMethod', 'kdtree', 'Distance', 'chebychev');
        fprintf('N');
        numChi2 = 0; % number of regional maxima merged into existing ROIs because they passed a Chi^2 test
        for i = 1:length(regmax)
            if assignment(i) == 0
                % Get region within watershed that overlaps other ROIs
                rng = ROIRng([xRegmax(i);yRegmax(i);zRegmax(i)]);
                region = false(size(watersheds));
                for j = 1:length(allNeighbors{i})
                    region(rng{:}) = region(rng{:}) | ROIShapes(:,:,:,allNeighbors{i}(j));
                end
                region = region & watersheds == i;

                % See if residual passes Chi^2 test
                if 1 - chi2cdf(norm(residual(region))^2,nnz(region)) > params.pval
                    % Assign regional maximum to ROI with greatest power
                    pow = zeros(length(allNeighbors{i}),1);
                    for ii = 1:length(allNeighbors{i})
                        rng = ROIRng(ROIOffset(:,allNeighbors{i}(ii)));
                        pow(ii) = norm(rates(allNeighbors{i}(ii))*vec(ROIShapes(region(rng{:}),allNeighbors{i}(ii))));
                    end
                    [~,jj] = max(pow);
                    assignment(i) = allNeighbors{i}(jj);
                    numChi2 = numChi2 + 1;
                end
            end
        end

        for i = 1:length(regmax)
            if assignment(i) == 0 % create new ROI
                if numROI > params.maxROI, pause; end

                [ROICenter(1,numROI), ROICenter(2,numROI), ROICenter(3,numROI)] = ind2sub(params.sz,regmax(i));
                ROIOffset(:,numROI) = int32(ROICenter(:,numROI));
                rng = ROIRng(ROIOffset(:,numROI));
                ROIShapes(:,:,:,numROI) = residual(rng{:}) .* (watersheds(rng{:})==i) / intensity(i);
                ROIPrec(:,:,:,numROI) = intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i));
                ROIPower(numROI) = intensity(i)^2;
                numROI = numROI + 1;
            else % merge ROI
                j = assignment(i);
                rng = ROIRng(ROIOffset(:,j));
                ROIShape(:,:,:,j) = (ROIPrec(:,:,:,j) .* ROIShape(:,:,:,j) + ...
                    intensity(i)*(watershed(rng{:})==i)/(params.var + params.varSlope*intensity(i)) .* data(rng{:})) ./ ... % this line right here might be why sometimes we get multiple ROIs mixed together. And why aren't we updating all ROIs?
                (ROIPrec(:,:,:,j) + intensity(i)^2*(watershed(rng{:})==i)/(params.var + params.varSlope*intensity(i)));
                ROICenter(:,j) = (ROIPower(j) * ROICenter(:,j) + intensity(i).^2 * [xRegmax(i); yRegmax(i), zRegmax(i)])...
                    /(ROIPower(j) + intensity(i).^2);
                ROIPower(j) = ROIPower(j) + intensity(i)^2;
                ROIPrec(:,:,:,j) = ROIPrec(:,:,:,j) + intensity(i)^2*(watershed(rng{:})==i)/(params.var + params.varSlope*intensity(i));
            end
        end
        ROIShape(isnan(ROIShape)) = 0;
    else
        numROI = numROI + length(regmax);
        for i = 1:numROI
            [ROICenter(1,i), ROICenter(2,i), ROICenter(3,i)] = ind2sub(params.sz,regmax(i));
            ROIOffset(:,i) = int32(ROICenter(:,i));
            rng = ROIRng(ROIOffset(:,i));
            ROIShapes(:,:,:,i) = data(rng{:}) .* (watersheds(rng{:})==i) / intensity(i);
            ROIPrec(:,:,:,i) = intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i));
        end
        ROIPower(1:numROI) = intensity.^2;
        numNeighbors = 0;
        numChi2 = 0;
    end
    fprintf('\t Found %d regional maxima: %d neighbors, %d pass Chi^2, %d new, %d total ROI, %gs\n', length(regmax), numNeighbors, numChi2, length(regmax)-numNeighbors-numChi2, numROI, toc);
end
