function [ tifinfo ] = get_cohen_info( tifdir )
%GET_COHEN_INFO Get the info about a tif stack, including:
%   - flist:    a list of filenames
%   - imsize:   dimensions of each frame
tifdir
flist = dir(fullfile(tifdir, '*.tif'));
nf = length(flist);
fnum = zeros(nf, 1);

% Sort the tifs by the index in their name 
% eg. fname-index.tif
for j = 1:nf;
    indx1 = strfind(flist(j).name, '-') + 1;
    indx2 = strfind(flist(j).name, '.tif') - 1;
    fnum(j) = str2num(flist(j).name(indx1:indx2));
end;
[~, findx] = sort(fnum);
flist = flist(findx);

% Get the size of the images from the first tif
tmp = imread(fullfile(tifdir, flist(1).name));
imsize = size(tmp);
nframes = length(flist);

% Package into struct
tifinfo.flist = flist;
tifinfo.imsize = imsize;
tifinfo.nframes = nframes;

end

