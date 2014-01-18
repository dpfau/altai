params.sig = 5;
params.dz = 12.5;
params.thresh = 0.012;

for t = 1
    data = loadframe(t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('S');
    clear gpuData; % save space on the GPU
    gpuRegmax = myregionalmax(gpuDataBlur-params.thresh);
    fprintf('R');
    gpuWatershed = arrayfun(@(seed) mywatershed(gpuDataBlur-params.thresh, seed), find(gpuRegmax), 'UniformOutput', 0);
    fprintf('W');
end
