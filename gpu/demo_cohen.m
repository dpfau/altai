clear params
addpath(fullfile('..','util'));
% tifdir = fullfile('/','Users','scott','Projects','cohen','data','coke_can_08_01');
tifdir = fullfile('C:','Users','scott','Projects','cohen','data','Take 1_Arch_2ms_exp', 'tifs');
tifdir = fullfile('/', 'group','hips','scott',...
                  'data','cohen','2013_08_01_coke_can_wt','tifs');
tifinfo = get_cohen_info(tifdir);

params.sig = [5,5];
params.dz = 1;
params.thresh = 1200;
params.roiSz = [60,60]; % Max size of an ROI

% Get 3-dim image size (x,y,1)
params.sz = tifinfo.imsize;
params.minDist = 3; % minimum distance below which two ROI are considered the same
params.autoVar = false; % automatically compute baseline variance and linear dependence on firing rate
params.var = 10;
params.varSlope = 0;
params.maxROI = 1e5;
params.pval = 1e-14; % very strict.
params.watershedFlag = false;
if exist('tRng','var')
    params.tRng = tRng;
else
    params.tRng = [1:10]; 
end

% params.loadframe = @(t) double(imread(fullfile(tifdir, tifinfo.flist(t).name)));
params.loadframe = @(t) load_cohen(t, tifdir, tifinfo);

% params.minSize = 30; % minimum number of pixels to consider in a watershed
% params.maxSize = prod(params.roiSz); % maximum number of pixels to consider in a watershed (set extremely high at the moment)

ROI = detectROIs(params);
save ROI_results ROI
