data_flag = 1;
sz = [400,150,100];
if ~exist('data','var')
    if data_flag == 0
        sz = [400,150,100];
        load ~/Dropbox/From' Pfau to Gabitto'/test_dat_s30_t1-100_smallcrop_xt.txt
        data = padarray(reshape(test_dat_s30_t1_100_smallcrop_xt,sz),[20,20]);
        clear test_dat_s30_t1_100_smallcrop_xt
    else
        sz = [150,250,10,100];
        load ~/Dropbox/From' Pfau to Gabitto'/test_dat_s25-34_t1-100_smallercrop_xt.txt
        data = reshape(test_dat_s25_34_t1_100_smallercrop_xt,sz);
        clear test_dat_s25_34_t1_100_smallercrop_xt;
        if data_flag == 1
            data = data(:,:,7,:);
        elseif data_flag == 2
            data = data(:,:,5,:);
        end
        data = padarray(data,[20,20]);
    end
    sz = size(data);
end

params.var_slope = .0008; % The slope of the variance-vs-mean relationship
params.var_offset = .0055; % The variance at 0 mean
params.blobsize = 6;
params.threshold = 1.2;
params.regmax_dist = 5;
params.min_size = 50; % minimum number of pixels in an ROI
params.pval = 1e-14; % if it's good enough for particle physics, it's not good enough for me. even stricter!
params.optose = true; % display ROIs at each step 
params.pause = false; % stop at each frame to show progress?

ROI = extract_roi_morph(data,params);