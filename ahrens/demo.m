params.blobsize = 5;
params.thresh = 0.012;
params.dz = 12.5;

for t = 1:1000
    fprintf('Loading frame %d...',t);
    data = loadframe(t);
    
    fprintf('Smoothing...');
    datablur = imgaussian3(data,params.blobsize,params.dz);
    
    fprintf('Regional max...')
    regmax = datablur .* imregionalmax(abs(datablur).*(abs(datablur)>params.thresh));
    
    fprintf('Watershedding...')
    wtrshed = uint16(watershed(-abs(datablur).*(abs(datablur)>params.thresh))).*uint16(abs(datablur)>params.thresh);
    
    fprintf('\n');
end