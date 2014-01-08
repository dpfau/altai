% sz = [400,150,100];
% if ~exist('data','var')
%     load ../../Ahrens/test_dat_s30_t1-100_smallcrop_xt.txt
%     data = padarray(reshape(test_dat_s30_t1_100_smallcrop_xt,sz),[20,20]);
%     clear test_dat_s30_t1_100_smallcrop_xt
%     sz = size(data);
% end

if ~exist('data','var')
    load ../../Ahrens/test_dat_s25-34_t1-100_smallercrop_xt.txt
    data = reshape(test_dat_s25_34_t1_100_smallercrop_xt,[150,250,10,100]);
    clear test_dat_s25_34_t1_100_smallercrop_xt;
    data = data(:,:,7,:);
    data = padarray(data,[20,20]);
    sz = size(data);
end

params.sig = 0.07;
params.mode = 'hermite-gaussian';
params.window = 41;
params.order = 8;
params.blobsize = 6;
params.threshold = 1.2;
params.optose = true; % display ROIs at each step 
params.pause = false; % stop at each frame to show progress?

ROI = extract_roi(data,params);