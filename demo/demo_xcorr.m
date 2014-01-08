sz = [400,150,100];
load ../../Ahrens/test_dat_s30_t1-100_smallcrop_xt.txt
data = reshape(test_dat_s30_t1_100_smallcrop_xt,sz);
clear test_dat_s30_t1_100_smallcrop_xt

filt = conj(fft2(make_blob(sz(1),sz(2),10)));
localmax = zeros(sz);
thresh = 10;
for i = 1:sz(3)
    smoothed_data = ifft2(fft2(data(:,:,i)).*filt);
    localmax(:,:,i) = abs(smoothed_data).*(abs(smoothed_data) > thresh & imregionalmax(abs(smoothed_data)));
end

% remove regional maxima on the border
localmax(1,:,:) = 0;
localmax(sz(1),:,:) = 0;
localmax(:,1,:) = 0;
localmax(:,sz(2),:) = 0;

% ROI = struct(); % struct array containing location and shape of ROI
% r = 10;
% maxIter = 10000;
% % Find shapes and locations of ROIs
% for t = 1:maxIter
%     [x,y] = ind2sub(sz(1:2),find(localmax(:,:,mod(t-1,100)+1)));
%     for i = 1:length(x)
%         k = 0; % index of this ROI in the list...if found. Note: this won't scale well! Something like KD-tree is way more efficient way to search.
%         for j = 1:length(ROI)
%             if (ROI(j).x-x)^2 + (ROI(j).y-y)^2 < r^2
%                 if k == 0
%                     k = j;
%                 else
%                     % conflict! This is an edge case, figure out what to do later.
%                 end
%             end
%             if k == 0
%                 k = length(ROI);
%                 ROI(k).x = x;
%                 ROI(k).y = y;
%             end
%         end
%     end
% end

% Frame-by-frame, fit the most likely shape/position/intensity for each ROI
eta = 0.07; % standard deviation of pixel noise
sigma = 1; % standard deviation of each shape from the template
xi = 10; % standard deviation of regmax position