function ROI = extract_roi_morph(data,params)
% A morphological, rather than statistical, approach to identifying ROIs

sz = size(data);
ROI = [];
params.T = sz(end);
for i = 1:sz(end)
    fprintf('\nAnalyzing frame %d of %d',i,sz(end));
    ROI = update_roi_morph(data(:,:,i),ROI,params,i);
end