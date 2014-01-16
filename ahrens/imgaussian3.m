function Y = imgaussian3(X,sig,dz)
% Gaussian blur of an image, where the kernel width in the z direction is
% significantly smaller (so much so that the kernel width is trivially
% small)
%
% Input:
%   X - 3D data to be blurred
%   sig - kernel width in the x and y directions
%   dz - ratio between the kernel width in the x/y direction and z
%       direction (so width in z direction is sig/dz)

X_ = zeros(size(X));
n = size(X,3);
for i = 1:n
    X_(:,:,i) = imgaussian(X(:,:,i),sig);
end

Y = zeros(size(X));
lim = floor(sqrt(-(sig/dz)^2*log(1e-6)));
for i = -lim:lim
    Y(:,:,max(0,i)+1:min(0,i)+n) = Y(:,:,max(0,i)+1:min(0,i)+n) + ...
        exp(-i^2/(sig/dz)^2)/sqrt(2*pi)/(sig/dz)*X_(:,:,max(0,-i)+1:min(0,-i)+n);
end