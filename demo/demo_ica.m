load ~/Documents/Research/Janelia/Ahrens/test_dat_s30_t1-100_smallcrop_xt.txt
data = reshape(test_dat_s30_t1_100_smallcrop_xt,400,150,100);
wtrsheds = zeros(size(data));
regmaxs = zeros(size(data));
for i = 1:100
    [regmax,wtrshed] = get_regmax(data(:,:,1),6,1.2);
    regmaxs(:,:,i) = regmax;
    wtrsheds(:,:,i) = wtrshed;
end

parts = zeros(400,150,nnz(regmaxs));
idx = find(regmaxs);
for i = 1:length(idx);
    [x,y,z] = ind2sub([400,150,100],idx(i));
    wtrshed_id = wtrsheds(x,y,z);
    parts(:,:,i) = data(:,:,z).*(wtrsheds(:,:,z)==wtrshed_id);
end
