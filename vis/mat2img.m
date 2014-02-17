function img = mat2img(X)
% Makes an image where each slice of X is given a random color
assert(ndims(X)==3)
[p,q,n] = size(X);
cols = hsv(n);
cols = cols(randperm(n),:);
img = zeros(p,q,3);
X = reshape(X,p*q,n);
for i = 1:n
    img = img + reshape(X(:,i)*cols(i,:),p,q,3);
end
img = mat2gray(img);