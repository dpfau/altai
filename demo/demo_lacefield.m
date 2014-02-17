if ~exist('data','var')
    data = zeros(135,137,3972);
    for i = 1:3972
        data(:,:,i) = imread('/Users/pfau/Dropbox/Pfau/2013-07-25-010.tif','Index',i);
    end
    sz = size(data);
end

params.var_slope = 0; % The slope of the variance-vs-mean relationship
params.var_offset = 6.44e5; % The variance at 0 mean
params.blobsize = 2;
params.threshold = 2e4;
params.regmax_dist = 5;
params.min_size = 50; % minimum number of pixels in an ROI
params.pval = 1e-14; % if it's good enough for particle physics, it's not good enough for me. even stricter!
params.optose = true; % display ROIs at each step 
params.pause = false; % stop at each frame to show progress?

ROI = extract_roi_morph(data,params);