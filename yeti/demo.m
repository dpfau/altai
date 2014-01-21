params.sig = 5;
params.dz = 12.5;
params.thresh = 0.011;
params.sz = [1472,2048,41];

wtr = zeros(params.sz,'int32'); % initialize up front for speed
for t = 1
    data = loadframe(t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('S');
    clear gpuData; % save space on the GPU
    gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
    fprintf('R');
    fastwatershed(gather(gpuDataBlur-params.thresh), wtr, gather(gpuRegmax));
    fprintf('W');
end
