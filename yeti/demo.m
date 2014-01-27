params.sig = 5;
params.dz = 12.5;
params.thresh = 0.011;
params.roiSz = [40,40,3]; % Max size of an ROI
params.sz = [1472,2048,41] + floor(params.roiSz/2); % don't forget the padding!

numROI = int32(0);
ROIShapes = zeros([params.roiSz,1e5]); % Initialize the whole sparse array. Expanding as we go is slow and dumb.
ROIPrecs  = zeros([params.roiSz,1e5]); % the precision of each pixel in the ROI.
ROIOffset = zeros(3,1e5,'int32'); % Location of the ROIs. Fixed from the start, integer precision
ROICenter = zeros(3,1e5); % Slightly different from ROIOffset. Double precision, updated online, used to decide if two ROI are close enough to merge.
for t = 1
	watersheds = zeros(params.sz,'int32');
    data = padarray(loadframe(t),floor(params.roiSz/2)); % Prevents ROIs from going over the edge of an image.
    fprintf('%d',t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('.');
    clear gpuData; % save space on the GPU
    gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
    regmax = gather(gpuRegmax);
    fprintf('.');
    fastwatershed(gather(gpuDataBlur-params.thresh), watersheds, regmax);
    fprintf('.');
    if numROI > 0
        residual = double(data); 
        rates = ratesfromframe(residual, ROIShapes, ROIOffset, numROI);
        getResidual(data,residual,ROIShapes,ROIOffset,rates);
        % Again, this part could be a bottleneck when the number of ROIs is large, and it may make more sense to use a KD-tree like datastructure to store the ROI locations
    else
        numROI = numROI + max(watersheds(:));
        for i = 1:numROI
            [ROICenter(1,i), ROICenter(2,i), ROICenter(3,i)] = ind2sub(params.sz,regmax(i));
            ROIOffset(:,i) = int32(ROICenter(:,i)-floor(params.roiSz/2)');
            rng = arrayfun(@(x,y)x+(1:int32(y)),ROIOffset(:,i),params.roiSz','UniformOutput',0);
            ROIShapes(:,:,:,i) = data(rng{:}).*(watersheds(rng{:})==regmax(i));
        end
    end
    fprintf('.');
    fprintf('\n');
end
