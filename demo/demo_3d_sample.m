clear params
params.sig = [5,5,5/12.5];
params.thresh = 0.03;
params.roiSz = [40,40,3]; % Max size of an ROI
params.sz = [150,250,12];
params.minDist = 3; % minimum distance below which two ROI are considered the same
params.autoVar = true; % automatically compute baseline variance and linear dependence on firing rate
params.maxROI = 1e5;
params.pval = 1e-14; % very strict.
params.tRng = 1:100;
params.watershedFlag = false;

load('/Users/pfau/Dropbox/Ahrens Lab Data Sample/3D movie sample/test_dat_s25-34_t1-100_smallercrop_xt.txt')
data = padarray(reshape(test_dat_s25_34_t1_100_smallercrop_xt,150,250,10,100),[0,0,1,0]);
clear test_dat_s25_34_t1_100_smallercrop_xt
params.loadframe = @(t) data(:,:,:,t);

% params.minSize = 30; % minimum number of pixels to consider in a watershed
% params.maxSize = prod(params.roiSz); % maximum number of pixels to consider in a watershed (set extremely high at the moment)

ROI = detectROIs(params);
save ROI_results_3d_sample ROI
