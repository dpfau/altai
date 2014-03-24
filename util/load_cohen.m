function [ frame ] = load_cohen( t, tifdir, tifinfo )
%LOAD_COHEN Load a tiff stack stored as a directory full of tifs named
% according to:
%   fname-index.tif
% where 'index' is an integer id for each image.

frame = double(imread(fullfile(tifdir, tifinfo.flist(t).name)));

% Subtract off the mean image
frame = frame - tifinfo.meanimg;

% Call ZCA
frame = tifinfo.zca(frame);

end

