function C = ROICov(ROIShapes,ROIOffset)
% Fast computation of covariance between ROI shapes

roiSz = size(ROIShapes);
patchrng = arrayfun(@(x)1:x,roiSz(1:end-1),'UniformOutput',0);

C = sparse(roiSz(end),roiSz(end));
warning('off','stats:KDTreeSearcher:rangesearch:DataConversion'); % don't need to hear about my conversions
neighbors = rangesearch(double(ROIOffset)'/diag(roiSz(1:end-1)), ...
                        double(ROIOffset)'/diag(roiSz(1:end-1)), 1, ...
                        'NSMethod', 'kdtree', 'Distance', 'chebychev');
for i = 1:roiSz(end)
    for j = 1:length(neighbors{i})
        k = neighbors{i}(j);
        C(i,k) = ROIDot(ROIShapes(patchrng{:},i),ROIShapes(patchrng{:},k),ROIOffset(:,i),ROIOffset(:,k));
    end
end

function z = ROIDot(x,y,px,py)

roiSz = size(x);
ROIRng = @(p1,p2) arrayfun(@(x,y,z) max(0,y-x)+1:min(0,y-x)+z,p1,p2,roiSz','UniformOutput',0);
vec = @(x)x(:);
rng1 = ROIRng(px,py);
rng2 = ROIRng(py,px);
z = vec(x(rng1{:}))'*vec(y(rng2{:}));