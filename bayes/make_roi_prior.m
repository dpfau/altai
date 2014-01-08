function [mu,Sig] = make_roi_prior(order)

% Prior mean and variance for each basis element. Set by hand based on empirical fiddling with actual ROI shapes.
[xx,yy] = meshgrid(0:order,0:order);
mu  = zeros(order+1);
Sig = .05./(xx+yy);
mu(1,1) = 1;
mu(1,2) = -.1;
mu(2,1) = -.1;
Sig(1,1) = 0.0002;
Sig(1,2) = 0.005;
Sig(2,1) = 0.005;

mu = mu(:);
Sig = diag(Sig(:));