function X = make_blob(h,w,r,kind)
% makes a rectangular image of size [h,w] with a blob of radius r in the
% upper left corner (and wrapping around the edge)

[y,x] = meshgrid(floor(-w/2):floor(w/2 - 1),floor(-h/2):floor(h/2-1));

if nargin < 4
    kind = 'gaussian';
end

switch kind
    case 'gaussian'
        X = exp(-(fftshift(x).^2 + fftshift(y).^2)/r^2)/sqrt(pi/2)/r; % normalized Gaussian
    case 'circular'
        X = sqrt(1 - (fftshift(x).^2 + fftshift(y).^2)/r^2);
        X(imag(X)~=0) = 0;
end