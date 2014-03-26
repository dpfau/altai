clear params
params.sig = 5;
params.dz = 12.5;
params.thresh = 0.006;
params.roiSz = [40,40,3]; % Max size of an ROI
params.sz = [1024,2048,43];
params.minDist = 3; % minimum distance below which two ROI are considered the same
params.autoVar = true; % automatically compute baseline variance and linear dependence on firing rate
params.maxROI = 1e5;
params.pval = 1e-14; % very strict.
if exist('tRng','var')
    params.tRng = tRng;
else
    params.tRng = [25:1000,1:24]; 
end

if gpuDeviceCount
    params.loadframe = @(t) padarray(loadframe(t),[0,0,1]);
else
    params.loadframe = @(t) padarray(loadframe(t,'/Users/pfau/Documents/Research/ROI/Light Sheet Data/binary_data'),[0,0,1]);
end
% params.minSize = 30; % minimum number of pixels to consider in a watershed
% params.maxSize = prod(params.roiSz); % maximum number of pixels to consider in a watershed (set extremely high at the moment)

ROI = detectROIs(params);
save ROI_results ROI
