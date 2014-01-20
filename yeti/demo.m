params.sig = 5;
params.dz = 12.5;
params.thresh = 0.012;

watershedKernel = parallel.gpu.CUDAKernel('watershedKernel.ptx','watershedKernel.cu');
watershedKernel.ThreadBlockSize = [512,1,1];
for t = 1
    data = loadframe(t);
    gpuData = gpuArray(data);
    gpuDataBlur = blur(gpuData, [params.sig, params.sig, params.sig/params.dz]);
    fprintf('S');
    clear gpuData; % save space on the GPU
    gpuRegmax = myregionalmax(gpuDataBlur-params.thresh);
    fprintf('R');
    gpuWatershed = feval(watershedKern, gpuDataBlur-params.thresh, uint32(find(regmax)));
    fprintf('W');
end
