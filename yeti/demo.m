params.sig = 5;
params.dz = 12.5;
params.thresh = 0.011;
params.sz = [1472,2048,41];
params.roiSz = [40,40,3]; % Max size of an ROI

numROI = int32(0);
ROIShapes = zeros([params.roiSz,1e5]); % Initialize the whole sparse array. Expanding as we go is slow and dumb.
ROIPrecs  = zeros([params.roiSz,1e5]); % the precision of each pixel in the ROI.
ROILocs   = zeros(3,1e5,'int32'); % location of the ROIs.
for t = 1
	watersheds = zeros(params.sz,'int32');
    data = loadframe(t);
    fprintf('%d',t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('.');
    clear gpuData; % save space on the GPU
    gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
    fprintf('.');
    fastwatershed(gather(gpuDataBlur-params.thresh), watersheds, gather(gpuRegmax));
    fprintf('.');
    if numROI > 0
        residual = double(data); 
        rates = ratesfromframe(residual, ROIShapes, ROILocs, numROI);
        getResidual(data,residual,ROIShapes,ROILocs,rates);
    else
        numROI = numROI + max(watersheds(:));
        for i = 1:numROI
            [ROILocs(1,i), ROILocs(2,i), ROILocs(3,i)] = ind2sub(params.sz,gpuRegmax(i))-floor(params.roiSz/2);
        end
    end
    fprintf('.');
    fprintf('\n');
end
