function intensity = get_intensity(frame,ROI,basis,intensity,fixed)

[m,n,q] = size(basis);
basis = reshape(basis,m*n,q);
for i = find(fixed)'
	xRng = ROI(i).pos(1)-floor(m/2)+(1:m);
    yRng = ROI(i).pos(2)-floor(n/2)+(1:n);
    frame(xRng,yRng) = frame(xRng,yRng) - intensity(i)*reshape(basis*ROI(i).mu,m,n);
end

basisFlat = [];
for i = find(~fixed)'
	xRng = ROI(i).pos(1)-floor(m/2)+(1:m);
    yRng = ROI(i).pos(2)-floor(n/2)+(1:n);
    roiShape = zeros(size(frame));
    roiShape(xRng,yRng) = reshape(basis*ROI(i).mu,m,n);
    basisFlat = [basisFlat, roiShape(:)];
end

intensity(~fixed) = pinv(basisFlat)*frame(:);