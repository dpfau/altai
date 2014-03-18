function [rates,residual] = ratesFromFrame(data,ROIShapes,ROIOffset,numROI)

roiSz = size(ROIShapes);
patchrng = arrayfun(@(x)1:x,roiSz(1:end-1),'UniformOutput',0);
% if numROI > 1e3 % While LSQR sounds nice, in practice it seems to be stupidly slow (at least the Matlab implementation)
%     [rates,flag] = lsqr(@(x,mode)localMultiply(x,mode,ROIShapes,ROIOffset,numROI,size(data),roiSz(1:3)'),data(:),1e-6,250);
% else
    C = ROICov(ROIShapes(patchrng{:},1:numROI),ROIOffset(:,1:numROI));
    warning('off','MATLAB:singularMatrix');
    rates = C\localMultiply(data(:),'transp',ROIShapes,ROIOffset,numROI,size(data),roiSz(1:end-1)', patchrng);
% end
if nargout == 2
    residual = data - reshape(localMultiply(rates,'notransp',ROIShapes,ROIOffset,numROI,size(data),roiSz(1:end-1)', patchrng),size(data));
end
end

function y = localMultiply(x,mode,ROIShapes,ROIOffset,numROI,sz,roiSz, patchrng)
vec = @(x)x(:);
ROIRng = @(x) arrayfun(@(x,y)x-int32(floor(y/2))+(0:int32(y)-1),x,roiSz,'UniformOutput',0);
if strcmpi(mode,'notransp')
    y = zeros(sz);
    for i = 1:numROI
        rng = ROIRng(ROIOffset(:,i));
        y(rng{:}) = y(rng{:}) + x(i)*ROIShapes(patchrng{:},i);
    end
    y = y(:);
elseif strcmpi(mode,'transp')
    x = reshape(x,sz);
    y = zeros(numROI,1);
    for i = 1:numROI
        rng = ROIRng(ROIOffset(:,i));
        y(i) = vec(x(rng{:}))'*vec(ROIShapes(patchrng{:},i));
    end
end
end