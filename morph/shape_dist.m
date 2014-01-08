function distmat = shape_dist(ROI)

distmat = zeros(length(ROI));
for i = 1:length(ROI)
    for j = i+1:length(ROI)
        distmat(i,j) = ROI(i).shape(:)'*ROI(j).shape(:)/norm(ROI(i).shape(:))/norm(ROI(j).shape(:));
    end
end
distmat = distmat + distmat';
distmat(1:length(ROI)+1:end) = 1;