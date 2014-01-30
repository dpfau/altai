
function data = loadframe(t,dataPath,verbose)

if nargin < 2 || isempty(dataPath)
    dataPath = '/vega/stats/users/dbp2112/ahrens/data';
end
if nargin < 3, verbose = false; end

data = zeros(1472,2048,41,'uint16');
for i = 1:41
    fid = fopen(sprintf('%s/dff_aligned_T%d_slice%d',dataPath,t,i),'r');
    slice = fread(fid,'uint16');
    fclose(fid);
    data(:,:,i) = reshape(slice,1472,2048);
    data = data(225:1248,:,:); % The animal is entirely within these bounds, so we can just throw out about 1/3rd of the data.
    if verbose, fprintf('.'); end
end
data = (single(data) - 15040)/2^16; % empirically, this seems to be where zero is
if verbose, fprintf('\n'); end
