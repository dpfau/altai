function [ tifinfo ] = get_cohen_info( tifdir, framestart, framestop )
%GET_COHEN_INFO Get the info about a tif stack, including:
%   - flist:    a list of filenames
%   - imsize:   dimensions of each frame
fprintf('Enumerating files in %s\n', tifdir);
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

assert(length(flist) > framestart);
assert(length(flist) >= framestop);

flist = flist(framestart:framestop);

flist(1)
flist(end)

% Get the size of the images from the first tif
tmp = imread(fullfile(tifdir, flist(1).name));
imsize = size(tmp);
nframes = length(flist);

% Compute the mean image
if exist(fullfile(tifdir, 'preprocessing.mat'), 'file') > 0
    load(fullfile(tifdir, 'preprocessing.mat'));
else
    fprintf('Computing the mean image and SVD\n');
    downsample = 1;
    frameinds = 1:downsample:nframes;
    nnframes = length(frameinds);
    imgs = zeros(horzcat(imsize,nnframes));
    for i = 1:nnframes
        t = frameinds(i);
        if mod(i,10) == 0
            fprintf('%d/%d\n', i, nnframes);
        end
        imgs(:,:,i) = double(imread(fullfile(tifdir, flist(t).name)));
    end
    fprintf('\n');

    % Compute the mean and the svd (for zca)
    meanimg = mean(imgs,3);
    imgs = reshape(bsxfun(@minus,imgs,meanimg),prod(imsize),nnframes);
    
    fprintf('Computing the SVD of the data\n');
    [u,s,~] = svd(imgs,0);

    % Save the meanimg and the svd for reuse
    save(fullfile(tifdir, 'preprocessing.mat'), 'meanimg', 'u', 's');
end
% Package into struct
tifinfo.flist = flist;
tifinfo.imsize = imsize;
tifinfo.nframes = nframes;
tifinfo.meanimg = meanimg;
tifinfo.zca = @(img) reshape(u*(diag(1./(diag(s)+1))*(u'*reshape(img,prod(imsize),1))),imsize);

keyboard


end

