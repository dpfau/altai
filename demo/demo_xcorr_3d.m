sz = [150,250,10,100];
if ~exist('data','var')
    load ../../Ahrens/test_dat_s25-34_t1-100_smallercrop_xt.txt;
    data = reshape(test_dat_s25_34_t1_100_smallercrop_xt,sz);
    clear test_dat_s25_34_t1_100_smallercrop_xt
end

regmax = get_regmax(data,6,2);

make_regmax_vid(data,regmax,5);
make_regmax_vid(data,regmax,7);

regmax_inds = find(regmax);

