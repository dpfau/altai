function data = blur(data,sig)
% Anisotropic Gaussian blur. Works fine on regular matrices, but really shines on GPU.
%
% David Pfau, 2014

if length(sig) == 1 % isotropic case
    sig = sig*ones(1,ndims(data));
end
assert(ndims(data)==length(sig));
n = ndims(data);

for i = 1:n
    x = size(data,i);
    idx = 1:n;
    idx(i) = 1;
    idx(1) = i;
    bump = gpuArray(exp(-fftshift(floor(-x/2):floor(x/2-1)).^2/sig(i)^2)'/sqrt(2*pi)/sig(i));
    data = bsxfun(@times,fft(data,[],i),permute(fft(bump),idx));
    data = ifft(data,[],i);
end 

data = real(data);