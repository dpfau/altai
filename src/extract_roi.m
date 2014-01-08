function ROI = extract_roi(data,params)

params.basis = make_basis(params.mode,params.order,params.window);
sz = size(data);
ROI = [];
for i = 1:sz(end)
    fprintf('\nAnalyzing frame %d of %d',i,sz(end));
    ROI = update_roi_2d(data(:,:,i),ROI,params);
end
plot_roi_shapes(ROI,params.basis,size(data));