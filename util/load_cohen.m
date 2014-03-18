function [ frame ] = load_cohen( t, tifdir, tifinfo )
%LOAD_COHEN Load a tiff stack stored as a directory full of tifs named
% according to:
%   fname-index.tif
% where 'index' is an integer id for each image.

frame = double(imread(fullfile(tifdir, tifinfo.flist(t).name)));

% HACK - double the frames to get a 3d data
% frame = repmat(frame, [1, 1, 2]);

end

