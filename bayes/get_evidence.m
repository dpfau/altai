function evidence = get_evidence(data,ROI,basis,intensity,sig,pos)

k = length(ROI);
[m,n,q] = size(basis);
xRng = pos(1)-floor(m/2)+(1:m);
yRng = pos(2)-floor(n/2)+(1:n);

data = data(xRng,yRng); % For the local approximation to the evidence, chop out everything but the area around the ROI
W = zeros(m,n,q,k);
mu = zeros(q,k);
Sig1 = zeros(q*k);
for i = 1:k
	xRng1 = max(1,ROI(i).pos(1)-pos(1)+1):min(m,ROI(i).pos(1)-pos(1)+m);
	yRng1 = max(1,ROI(i).pos(2)-pos(2)+1):min(n,ROI(i).pos(2)-pos(2)+n);
	xRng2 = max(1,pos(1)-ROI(i).pos(1)+1):min(m,pos(1)-ROI(i).pos(1)+m);
	yRng2 = max(1,pos(2)-ROI(i).pos(2)+1):min(n,pos(2)-ROI(i).pos(2)+n);
	for j = 1:q
		W(xRng1,yRng1,j,i) = intensity(i)*basis(xRng2,yRng2,j);
	end
	mu(:,i) = ROI(i).mu;
	Sig((i-1)*q+(1:q),(i-1)*q+(1:q)) = ROI(i).Sig;
end
W = reshape(W,m*n,q*k);
mu = mu(:);
data = data(:) - W*mu;
Sig = sig^2*eye(m*n) + W*Sig*W';

evidence = -1/2*data'/Sig*data - sum(log(diag(chol(Sig))));