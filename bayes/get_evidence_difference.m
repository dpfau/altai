function evidence = get_evidence_difference(data,ROI,basis,intensity,sig,idx)

k = length(ROI);
[m,n,q] = size(basis);
xRng = ROI(idx).pos(1)-floor(m/2)+(1:m);
yRng = ROI(idx).pos(2)-floor(n/2)+(1:n);

data = data(xRng,yRng); % For the local approximation to the evidence, chop out everything but the area around the ROI
W = zeros(m,n,q,k);
mu = zeros(q,k);
Sig1 = zeros(q*k);
for i = 1:k
	xRng1 = max(1,ROI(i).pos(1)-ROI(idx).pos(1)+1):min(m,ROI(i).pos(1)-ROI(idx).pos(1)+m);
	yRng1 = max(1,ROI(i).pos(2)-ROI(idx).pos(2)+1):min(n,ROI(i).pos(2)-ROI(idx).pos(2)+n);
	xRng2 = max(1,ROI(idx).pos(1)-ROI(i).pos(1)+1):min(m,ROI(idx).pos(1)-ROI(i).pos(1)+m);
	yRng2 = max(1,ROI(idx).pos(2)-ROI(i).pos(2)+1):min(n,ROI(idx).pos(2)-ROI(i).pos(2)+n);
	for j = 1:q
		W(xRng1,yRng1,j,i) = intensity(i)*basis(xRng2,yRng2,j);
	end
	mu(:,i) = ROI(i).mu;
	Sig((i-1)*q+(1:q),(i-1)*q+(1:q)) = ROI(i).Sig;
end
W = reshape(W,m*n,q*k);
mu = mu(:);
data = data(:) - W*mu;
Sig2 = sig^2*eye(m*n) + W*Sig*W';

evidence = zeros(k,1);
evidence(idx) = -1/2*data'/Sig2*data - sum(log(diag(chol(Sig2))));

% Loop over the other ROIs, computing the evidence when you 
for i = 1:k
	if i ~= idx
		data_ = data + W(:,(idx-1)*q+(1:q))*mu((idx-1)*q+(1:q)); % Add the ROI we're removing back to the data
		W_ = W;
		W_(:,(i-1)*q+(1:q)) = W_(:,(i-1)*q+(1:q))*intensity(idx)/intensity(i); % scale the intensity of this ROI
        data_ = data_ - W_(:,(i-1)*q+(1:q))*mu((i-1)*q+(1:q)); % And subtract the rescaled ROI off the data
		W_ = W_(:,[1:(idx-1)*q,idx*q+1:end]);
		Sig_ = sig^2*eye(m*n) + W_*Sig1([1:(idx-1)*q,idx*q+1:end],[1:(idx-1)*q,idx*q+1:end])*W_';
		evidence(i) = -1/2*data_'/Sig_*data_ - sum(log(diag(chol(Sig_))));
	end
end