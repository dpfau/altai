load frame
x = -5:0.25:5;
hg = @(n) hermite(n,x).*exp(-x.^2/2)/norm(hermite(n,x).*exp(-x.^2/2)); % Hermite-Gaussian functions

k = length(x);
max_order = 15;
frame_approx = zeros([size(frame),max_order]);
for i = 1:max_order % approximate the frame with Hermite-Gaussian polynomials up to order 15
    hgmat = []; % matrix of 2D Hermite-Gaussian functions
    for p = 0:i
        for q = 0:i
            hgmat = cat(3,hgmat,hg(p)'*hg(q));
        end
    end
    hgmat = reshape(hgmat,k^2,(i+1)^2);
    frame_approx(:,:,i) = reshape(hgmat/(hgmat'*hgmat)*hgmat'*vec(frame),k,k);
end
figure(1)
imagesc(frame); axis image

st = @(x,t) sign(x).*max(0,abs(x)-t); % Soft threshold, to clean up some junk
figure(2)
for i = 1:max_order
    subplot(2,max_order,i)
    imagesc(frame_approx(:,:,i)), axis image
    subplot(2,max_order,i+max_order)
    imagesc(st(frame_approx(:,:,i),.1)), axis image
end