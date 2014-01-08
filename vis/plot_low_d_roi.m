function plot_low_d_roi(ROI,basis,sz)

st = @(x,t) sign(x).*max(0,abs(x)-t); % soft threshold the ROI shape to clean it up a little
k = length(ROI);
cols = hsv(k);
cols = cols(randperm(k),:);
img = zeros(prod(sz(1:2)),3);
foo = zeros(sz(1:2));
[m,n,q] = size(basis);
for i = 1:k
    roi_shape = reshape(reshape(basis,m*n,q)*ROI(i).mu,m,n);
    offset = ROI(i).pos - floor([m,n]/2);
    foo(offset(1)+(1:m),offset(2)+(1:n)) = roi_shape;
    img = img + st(foo(:)/norm(foo),0.02)*cols(i,:);
    foo(offset(1)+(1:m),offset(2)+(1:n)) = 0;
end

figure(62)
clf
image(mat2gray(reshape(img,[sz(1:2),3])))
axis image
hold on
for i = 1:k
    scatter(ROI(i).pos(2),ROI(i).pos(1),100,'ro','LineWidth',2)
end