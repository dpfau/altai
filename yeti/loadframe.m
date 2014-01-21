
function data = loadframe(t,dataPath)

if nargin < 2
    dataPath = '/vega/stats/users/dbp2112/ahrens/data';
end
data = zeros(1472,2048,41,'uint16');
for i = 1:41
    fid = fopen(sprintf('%s/dff_aligned_T%d_slice%d',dataPath,t,i),'r');
    slice = fread(fid,'uint16');
    fclose(fid);
    data(:,:,i) = reshape(slice,1472,2048);
    fprintf('.');
end
data = (single(data) - 15040)/2^16; % empirically, this seems to be where zero is
fprintf('\n');
