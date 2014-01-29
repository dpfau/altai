params.sig = 5;
params.dz = 12.5;
params.thresh = 0.011;
params.roiSz = [40,40,3]; % Max size of an ROI
params.sz = [1472,2048,41] + 2*floor(params.roiSz/2); % don't forget the padding!
params.minDist = 3; % minimum distance below which two ROI are considered the same
params.var = 1; % baseline variance
params.varSlope = 0.001; % slope of the variance as a function of the intensity
params.maxROI = 1e5;

numROI = int32(0);
ROIShapes = zeros([params.roiSz,params.maxROI]); % Initialize the whole sparse array. Expanding as we go is slow and dumb.
ROIPrecs  = zeros([params.roiSz,params.maxROI]); % the precision of each pixel in the ROI.
ROIOffset = zeros(3,params.maxROI,'int32'); % Location of the ROIs. Fixed from the start, integer precision
ROICenter = zeros(3,params.maxROI); % Slightly different from ROIOffset. Double precision, updated online, used to decide if two ROI are close enough to merge.
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
        getResidual(data,residual,ROIShapes,ROIOffset,rates); fprintf('R');

        % Compute nearest neighbors, if regional maxima are close enough to ROI centers, merge them together
        [xRegmax, yRegmax, zRegmax] = ind2sub(params.sz,regmax);
        zRegmax = zRegmax * params.dz;
        warning('off','stats:KDTreeSearcher:knnsearch:DataConversion');
        [nearestNeighbors, nnDistance] = knnsearch((diag([1,1,params.dz])*ROICenter(:,1:numROI))', [xRegmax,yRegmax,zRegmax], 'K', 1);
        assignment(nnDistance < params.minDist) = nearestNeighbors(nnDistance < params.minDist);
        numNeighbors = nnz(nnDistance < params.minDist);

        % Get index of ROIs that overlap regional maxima
        warning('off','stats:KDTreeSearcher:rangesearch:DataConversion'); % don't need to hear about my conversions
        allNeighbors = rangesearch((diag([1,1,params.dz])*ROICenter(:,1:numROI))',[xRegmax,yRegmax,zRegmax], max(params.sz .* [1,1,params.dz]), 'NSMethod', 'kdtree', 'Distance', 'chebychev');
        fprintf('N');
        numResiduals = 0;
        for i = 1:length(regmax)
            if assignment(i) == 0

            end
        end

        for i = 1:length(regmax)
            if assignment(i) == 0 % create new ROI
                if numROI > params.maxROI, pause; end

                [ROICenter(1,numROI), ROICenter(2,numROI), ROICenter(3,numROI)] = ind2sub(params.sz,regmax(i));
                ROIOffset(:,numROI) = int32(ROICenter(:,numROI));
                rng = arrayfun(@(x,y)x-floor(int32(y)/2)+(0:int32(y)-1),ROIOffset(:,numROI),params.roiSz','UniformOutput',0);
                ROIShapes(:,:,:,numROI) = residual(rng{:}) .* (watersheds(rng{:})==i) / intensity(i);
                ROIPrec(:,:,:,numROI) = intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i));
                numROI = numROI + 1;
            else % merge ROI
            end
        end
    else
        numROI = numROI + length(regmax);
        for i = 1:numROI
            [ROICenter(1,i), ROICenter(2,i), ROICenter(3,i)] = ind2sub(params.sz,regmax(i));
            ROIOffset(:,i) = int32(ROICenter(:,i));
            rng = arrayfun(@(x,y)x-floor(int32(y)/2)+(0:int32(y)-1),ROIOffset(:,i),params.roiSz','UniformOutput',0);
            ROIShapes(:,:,:,i) = data(rng{:}) .* (watersheds(rng{:})==i) / intensity(i);
            ROIPrec(:,:,:,i) = intensity(i)^2 .* (watersheds(rng{:})==i) / (params.var + params.varSlope*intensity(i));
            numNeighbors = 0;
            numResiduals = 0;
        end
    end
    fprintf('\t Found %d regional maxima: %d neighbors, %d residuals, %d new, %d total ROI, %gs\n', length(regmax), numNeighbors, numResiduals, length(regmax)-numNeighbors-numResiduals, numROI, toc);
end
