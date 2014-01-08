figure(1)
clf
colormap gray
imagesc(frame)
axis image
hold on
scatter(y,x,50,'g.')
for i = 1:k
    text(y(i)+10,x(i),num2str(i),'Color','g');
end
hold off

resid = frame; % residual

st = @(x,t) sign(x).*max(0,abs(x)-t); % soft threshold the ROI shape to clean it up a little
k = length(ROI);
cols = hsv(k);
cols = cols(randperm(k),:);
img = zeros(numel(frame),3);
foo = zeros(size(frame));
[m,n,q] = size(params.basis);
for i = 1:k
    roi_shape = reshape(reshape(params.basis,m*n,q)*ROI(i).mu,m,n);
    offset = ROI(i).pos - floor([m,n]/2);
    foo(offset(1)+(1:m),offset(2)+(1:n)) = roi_shape;
    img = img + st(foo(:)/norm(foo),0.02)*cols(i,:);
    resid = resid - regmax(x(i),y(i))*foo;
    foo(offset(1)+(1:m),offset(2)+(1:n)) = 0;
end

figure(2)
clf
colormap gray
imagesc(resid)
axis image

figure(3)
clf
image(mat2gray(reshape(img,[size(frame),3])))
axis image

drawnow