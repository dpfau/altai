function [mu,alpha] = get_roi_posterior(frame,regmax,sig1,sig2,eta1,eta2)
% Get the mean ROI shape from a Gaussian prior and knowledge of the
% regularized maximum mean and intensity. Leaving out the covariance for
% now because it's too high-dimensional, will probably have to take some
% kind of approximation such as only keeping the diagonal (intuitively,
% there should really only be cross terms in the posterior covariance if
% two ROIs overlap)
%
% frame - one frame of data
% regmax - the location and intensity of regularized maxima, same size as
%   frame
% sig1 - the width of the gaussian bump in the mean ROI
% sig2 - the width of the gaussian bump in the ROI covariance (which is
%   diagonal)
% eta1 - the noise variance in the image
% eta2 - the width of the ROI shape variance (i.e. height of the Gaussian
%   bump)
%
% David Pfau, 2013

assert(~any(size(frame)~=size(regmax)));
[x,y] = ind2sub(size(regmax),find(regmax));
k = length(x);
mu = zeros([size(frame),k]);
Sig = mu;
SigMarg = eta1*ones(size(frame)); % diagonal of the marginal covariance of the frame
alpha = zeros(k,1);
for i = 1:k
    alpha(i) = regmax(x(i),y(i)); % intesity of regularized maximum
    mu(:,:,i)  =      circshift(make_blob(size(frame,1),size(frame,2),sig1),[x(i),y(i)]-1);
    Sig(:,:,i) = eta2*circshift(make_blob(size(frame,1),size(frame,2),sig2,'circular'),[x(i),y(i)]-1);
    frame = frame - alpha(i)*mu(:,:,i);
    SigMarg = SigMarg + alpha(i)^2*Sig(:,:,i);
end
frame = frame./SigMarg;
for i = 1:k
    mu(:,:,i) = mu(:,:,i) + alpha(i)*Sig(:,:,i).*frame;
end