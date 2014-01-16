function X = make_blob(sz,sig,kind)
% makes a rectangular image of size [h,w] with a blob of radius r in the
% upper left corner (and wrapping around the edge)

assert(length(sig)==length(sz)||length(sig)==1);
if length(sig) == 1
    sig = sig*ones(size(sz));
end

idx = arrayfun(@(x) floor(-x/2):floor(x/2-1), sz, 'UniformOutput', 0);
if length(sz)==2
    [y,x] = meshgrid(idx{[2,1]});
    Y = {x,y};
elseif length(sz)==3
    [y,x,z] = meshgrid(idx{[2,1,3]});
    Y = {x,y,z};
else
    error('The designers of MATLAB were clearly only capable of thinking in 3 dimensions, tops')
end

if nargin < 4
    kind = 'gaussian';
end

Z = arrayfun(@(i) fftshift(Y{i}).^2/sig(i)^2, 1:length(sz), 'UniformOutput', 0);
Z = apply(@plus,Z{:});
switch kind
    case 'gaussian'
        X = exp(-Z)/sqrt(pi/2*prod(sig)); % normalized Gaussian
    case 'circular'
        X = sqrt(1 - Z);
        X(imag(X)~=0) = 0;
end