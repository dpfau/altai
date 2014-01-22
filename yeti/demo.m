params.sig = 5;
params.dz = 12.5;
params.thresh = 0.011;
params.sz = [1472,2048,41];

watersheds = zeros(params.sz,'int32'); % initialize up front for speed

numROI = 0;
ROIShapes = sparse(prod(params.sz),1e5); % Initialize the whole sparse array. Expanding as we go is slow and dumb.
ROIPrecs  = sparse(prod(params.sz),1e5); % the precision of each pixel in the ROI.
for t = 1
    data = loadframe(t);
    fprintf('%d',t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('.');
    clear gpuData; % save space on the GPU
    gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
    fprintf('.');
    fastwatershed(gather(gpuDataBlur-params.thresh), watersheds, gather(gpuRegmax));
    fprintf('.');
    [rates, residual] = ratesfromframe(double(data), ROIShapes, int32(numROI));
    fprintf('.');
    fprintf('\n');
end
