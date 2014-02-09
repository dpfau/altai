function C = ROICov(ROIShapes,ROIOffset)
% Fast computation of covariance between ROI shapes

sz = size(ROIShapes);
C = sparse(sz(4),sz(4));
warning('off','stats:KDTreeSearcher:rangesearch:DataConversion'); % don't need to hear about my conversions
neighbors = rangesearch(double(ROIOffset)'/diag(sz(1:3)), double(ROIOffset)'/diag(sz(1:3)), 1, 'NSMethod', 'kdtree', 'Distance', 'chebychev');
for i = 1:sz(4)
    for j = 1:length(neighbors{i})
        k = neighbors{i}(j);
        C(i,k) = ROIDot(ROIShapes(:,:,:,i),ROIShapes(:,:,:,k),ROIOffset(:,i),ROIOffset(:,k));
    end
end

function z = ROIDot(x,y,px,py)

sz = size(x);
ROIRng = @(p1,p2) arrayfun(@(x,y,z) max(0,y-x)+1:min(0,y-x)+z,p1,p2,sz','UniformOutput',0);
vec = @(x)x(:);
rng1 = ROIRng(px,py);
rng2 = ROIRng(py,px);
z = vec(x(rng1{:}))'*vec(y(rng2{:}));