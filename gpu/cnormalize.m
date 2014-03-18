function [ output_args ] = cnormalize( input_args )
% Anisotropic Gaussian blur. Works fine on regular matrices, but really shines on GPU.
%
% David Pfau, 2014

if length(sig) == 1 % isotropic case
    sig = sig*ones(1,ndims(data));
end
assert(ndims(data)==length(sig));
n = ndims(data);

datablur = data;
for i = 1:n
    if sig(i)
        x = size(data,i);
        idx = 1:n;
        idx(i) = 1;
        idx(1) = i;
        bump = exp(-fftshift(floor(-x/2):floor(x/2-1)).^2/sig(i)^2)'/sqrt(2*pi)/sig(i);
        if isa(data,'gpuArray')
	       bump = gpuArray(bump);
	   end
        datablur = ifft(bsxfun(@times,fft(datablur,[],i),permute(fft(bump),idx)),[],i);
    end
end 

datablur = real(datablur);