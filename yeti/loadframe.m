
function data = loadframe(t)

data = zeros(1472,2048,41,'uint16');
for i = 1:41
    fid = fopen(sprintf('/vega/stats/users/dbp2112/ahrens/data/dff_aligned_T%d_slice%d',t,i),'r');
    slice = fread(fid,'uint16');
    fclose(fid);
    data(:,:,i) = reshape(slice,1472,2048);
    fprintf('.');
end
data = (single(data) - 15040)/2^16; % empirically, this seems to be where zero is
fprintf('\n');
