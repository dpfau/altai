function p = chi_sq_test(x1,x2,prec1,prec2)
% Chi squared test to see if two samples x1 and x2 from Gaussian
% distributions with known precisions prec1 and prec2 (in this case
% diagonal, for simplicity) are likely to come from a distribution with the
% same mean

x = (prec1.*x1 + prec2.*x2)./(prec1 + prec2); % ML estimate of the mean


% rescale error so it should be N(0,I) distributed under null hypothesis
y1 = sqrt(prec1).*(x1-x); 
y2 = sqrt(prec2).*(x2-x);
z = y1(prec1~=0)'*y1(prec1~=0) + y2(prec2~=0)'*y2(prec2~=0);
p = 1 - chi2cdf(z,nnz(prec1)+nnz(prec2));