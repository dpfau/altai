function posterior = update_posterior(data,prior,basis,intensity,sig)
% data - Should be self-explanatory. One frame of data used to update the
%   posterior
% prior - a struct array of parameters for each ROI. Each entry in the
%   array has the field "mu" and "Sig", which are the mean and covariance
%   parameters for that ROI. The posterior covariance may have some cross
%   terms between ROIs but we set those to zero. Technically this is a
%   variational approximation, albeit a trivial one.
% basis - a matrix of basis elements for the shape of a neuron in one ROI.
%   Could be Fourier modes, or Hermite-Gaussian modes, or any other
%   orthogonal basis
% position - the position of each ROI center
% intensity - the magnitude of the regional maximum at this point
% sig - the standard deviation of the pixel noise

[m,n,q] = size(basis);
[mm,nn] = size(data);
basisFlat = reshape(basis,m*n,q);
assert(length(prior)==length(intensity));
k = length(prior); % number of ROIs
posterior = [];
iiii = [];
jjjj = [];
for i = 1:k
    xRng = prior(i).pos(1)-floor(m/2)+(1:m);
    yRng = prior(i).pos(2)-floor(n/2)+(1:n);
    data(xRng,yRng) = data(xRng,yRng) - intensity(i)*reshape(basisFlat*prior(i).mu,m,n);
    
    [jj,ii] = meshgrid(yRng,xRng);
    iii = [];
    jjj = [];
    for j = 1:q
        iii = cat(3,iii,sub2ind([mm,nn],ii,jj));
        jjj = cat(3,jjj,(i-1)*q+j*ones(m,n));
    end
    iiii = cat(4,iiii,iii);
    jjjj = cat(4,jjjj,jjj);
end
basisFull = sparse(iiii(:),jjjj(:),vec(tprod(basis,[1,2,3],intensity,4)),mm*nn,k*q);

WTW = basisFull'*basisFull;
Sig = WTW/sig^2;
for i = 1:k
    Sig((i-1)*q+(1:q),(i-1)*q+(1:q)) = Sig((i-1)*q+(1:q),(i-1)*q+(1:q)) + inv(prior(i).Sig);
end
dmu = basisFull'*data(:)/sig^2 - WTW/Sig*basisFull'*data(:)/sig^4;
dSig = WTW/Sig*WTW;

for i = 1:k
    posterior(i).pos = prior(i).pos;
    posterior(i).mu  = prior(i).mu + prior(i).Sig*dmu((i-1)*q+(1:q));
    posterior(i).Sig = prior(i).Sig - prior(i).Sig*WTW((i-1)*q+(1:q),(i-1)*q+(1:q))*prior(i).Sig/sig^2 ...
        + prior(i).Sig*dSig((i-1)*q+(1:q),(i-1)*q+(1:q))*prior(i).Sig/sig^4;
    posterior(i).nFrames = prior(i).nFrames + 1;
    posterior(i).maxIntensity = max(prior(i).maxIntensity,abs(intensity(i)));
end