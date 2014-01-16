params.sig = 5;
params.dz = 12.5;
params.thresh = 0.012;

for t = 1
    data = loadframe(t);
    gpuData = gpuArray(data);
    gpuDataBlur = gpuBlur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('S');
    clear gpuData; % save space on the GPU
    gpuRegmax = myregionalmax(gpuDataBlur-params.thresh);
    fprintf('R');
    gpuWatershed = mywatershed(gpuDataBlur-params.thresh);
    fprintf('W');
end
