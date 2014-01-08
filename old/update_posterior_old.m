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
Sig = sig^2*speye(mm*nn);
for i = 1:k
    xRng = prior(i).pos(1)-floor(m/2)+(1:m);
    yRng = prior(i).pos(2)-floor(n/2)+(1:n);
    data(xRng,yRng) = data(xRng,yRng) - intensity(i)*reshape(basisFlat*prior(i).mu,m,n);
    
    basisFull{i} = [];
    [jj,ii] = meshgrid(yRng,xRng);
    for j = 1:q
        basisFull{i} = [basisFull{i}, reshape(sparse(ii,jj,basis(:,:,j),mm,nn),mm*nn,1)];
    end
    Sig = Sig + intensity(i)^2*basisFull{i}*sparse(prior(i).Sig)*basisFull{i}';
end

SigData = Sig\data(:);

for i = 1:k
    posterior(i).pos = prior(i).pos;
    posterior(i).mu  = prior(i).mu + intensity(i)*prior(i).Sig*basisFull{i}'*SigData;
    posterior(i).Sig = prior(i).Sig;% - intensity(i)^2*basisFull{i}'/Sig*basisFull{i};
end

% The implementation below uses the Woodbury lemma but is definitely buggy.
% As usual, I should have tried the straightforward way first.

% dataProj = zeros(k*q,1);
% for i = 1:k
%     dataProj((i-1)*q+(1:q)) = intensity(i)*basisFlat'*reshape(data(xRng,yRng),m*n,1);
% end
% 
% SigInv = sparse(k*q,k*q);
% for i = 1:k
%     for j = 1:k
%         if i == j
%         	SigInv((i-1)*q+(1:q),(i-1)*q+(1:q)) = prior(i).Sig + intensity(i)^2/sig^2*eye(q); % nice thing about the basis we're using is it's orthonormal
%         else
%             if ~any(abs(prior(i).pos-prior(j).pos) > [m,n])
%                 xRng1 = max(0,prior(j).pos(1)-prior(i).pos(1))+1:m+min(0,prior(j).pos(1)-prior(i).pos(1));
%                 yRng1 = max(0,prior(j).pos(2)-prior(i).pos(2))+1:n+min(0,prior(j).pos(2)-prior(i).pos(2));
%                 xRng2 = max(0,prior(i).pos(1)-prior(j).pos(1))+1:m+min(0,prior(i).pos(1)-prior(j).pos(1));
%                 yRng2 = max(0,prior(i).pos(2)-prior(j).pos(2))+1:n+min(0,prior(i).pos(2)-prior(j).pos(2));
%                 SigInv((i-1)*q+(1:q),(j-1)*q+(1:q)) = ...
%                     intensity(i)*intensity(j)/sig^2*...
%                     reshape(basis(xRng1,yRng1,:),numel(xRng1)*numel(yRng1),q)'*...
%                     reshape(basis(xRng2,yRng2,:),numel(xRng2)*numel(yRng2),q);
%             end
%         end
%     end
% end
% dataProj  = SigInv\dataProj;
% dataProj2 = zeros(m,n);
% for i = 1:k
%     xRng = prior(i).pos(1)-floor(m/2)+(1:m);
%     yRng = prior(i).pos(2)-floor(n/2)+(1:n);
%     dataProj2 = dataProj2 + intensity(i)*reshape(basisFlat*dataProj((i-1)*q+(1:q)),m,n);
% end
% 
% for i = 1:k
%     xRng = prior(i).pos(1)-floor(m/2)+(1:m);
%     yRng = prior(i).pos(2)-floor(n/2)+(1:n);
%     posterior(i).mu = prior(i).mu + ...
%         intensity(i)*basisFlat'*reshape(data(xRng,yRng),m*n,1)/sig^2 - ...
%         intensity(i)*basisFlat'*reshape(dataProj2,m*n,1)/sig^4;
%     posterior(i).pos = prior(i).pos;
% end