function [rates,residual] = ratesFromFrame(data,ROIShapes,ROIOffset,numROI,roiSz)

rates = lsqr(@(x,mode)localMultiply(x,mode,ROIShapes,ROIOffset,numROI,size(data)),data(:));
if nargout == 2
    residual = data - localMultiply(rates,'notransp',ROIShapes,ROIOffset,numROI,size(data),roiSz);
end

function y = localMultiply(x,mode,ROIShapes,ROIOffset,numROI,sz,roiSz)

vec = @(x)x(:);
ROIRng = @(x) arrayfun(@(x,y)x-int32(floor(y/2))+(0:int32(y)-1),x,roiSz','UniformOutput',0);
if strcmpi(mode,'notransp')
    y = zeros(sz);
    for i = 1:numROI
        rng = ROIRng(ROIOffset(:,i));
        y(rng{:}) = y(rng{:}) + x(i)*ROIShapes(:,:,:,i);
    end
    y = y(:);
elseif strcmip(mode,'transp')
    x = reshape(x,sz);
    y = zeros(numROI,1);
    for i = 1:numROI
        rng = ROIRng(ROIOffset(:,i));
        y(i) = vec(x(rng{:}))'*vec(ROIShapes(:,:,:,i));
    end
end