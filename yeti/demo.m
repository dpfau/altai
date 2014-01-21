params.sig = 5;
params.dz = 12.5;
params.thresh = 0.012;

watershedKernel = parallel.gpu.CUDAKernel('watershedKernel.ptx','watershedKernel.cu');
watershedKernel.ThreadBlockSize = [512,1,1];
for t = 101:150
    data = loadframe(t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('S');
    clear gpuData; % save space on the GPU
    gpuRegmax = int32(find(myregionalmax(gpuDataBlur-params.thresh)));
    fprintf('R');
    gpuWatershed = gpuArray.zeros(size(gpuDataBlur),'int32');
    feval(watershedKernel, gpuDataBlur-params.thresh, gpuWatershed, gpuRegmax, length(gpuRegmax), ndims(data), size(data));
    fprintf('W');
end
