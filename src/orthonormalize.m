function [Y,M] = orthonormalize(X)
% Constructs an orthonormal matrix Y with the same span as X, 
% and an upper triangular matrix M such that Y = X*M

[m,n] = size(X);
M = sparse(n,n);
Y = sparse(m,n);
for i = 1:n
	M(1:i-1,i) = Y(:,1:i-1)'*X(:,i);
	foo = X(:,i) - Y(:,1:i-1)*M(1:i-1,i);
	M(i,i) = norm(foo);
	Y(:,i) = foo/M(i,i);
    fprintf('%d\n',i);
end