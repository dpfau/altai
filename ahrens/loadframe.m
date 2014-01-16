function data = loadframe(t)

data = zeros(1472,2048,41,'uint16');
for i = 1:41
    fid = fopen(sprintf('/Users/pfau/Documents/Research/ROI/Light Sheet Data/binary_data/dff_aligned_T%d_slice%d',t,i),'r');
    slice = fread(fid,'uint16');
    fclose(fid);
    data(:,:,i) = reshape(slice,1472,2048);
end
data = (single(data) - 15040)/2^16; % empirically, this seems to be where zero is