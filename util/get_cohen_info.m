function [ tifinfo ] = get_cohen_info( tifdir )
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

% Get the size of the images from the first tif
tmp = imread(fullfile(tifdir, flist(1).name));
imsize = size(tmp);
nframes = length(flist);

% Compute the mean image
if exist(fullfile(tifdir, 'meanimg.mat'), 'file') > 0
    load(fullfile(tifdir, 'meanimg.mat'));
else
    fprintf('Computing the mean image\n');
    downsample = 10;
    nnframes = length(2:downsample:nframes);
    imgs = zeros(horzcat(imsize,nnframes));
    i = 0;
    for t = 2:downsample:nframes
        i = i+1;
        if mod(t,100) == 0
            fprintf('%d/%d\n', t, nframes);
        end
        imgs(:,:,i) = double(imread(fullfile(tifdir, flist(t).name)));
    end
    fprintf('\n');

    meanimg = mean(imgs,3);
    imgs = reshape(bsxfun(@minus,imgs,meanimg),prod(imsize),nnframes);
    [u,s,v] = svd(imgs,0);

    % Save the meanimg for reuse
    save(fullfile(tifdir, 'meanimg.mat'), 'meanimg');
end
% Package into struct
tifinfo.flist = flist;
tifinfo.imsize = imsize;
tifinfo.nframes = nframes;
tifinfo.meanimg = meanimg;
tifinfo.zca = @(img) reshape(u*(diag(1./diag(s))*(u'*reshape(img,prod(imsize),1))),imsize);

end

