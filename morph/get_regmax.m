function [regmax,wtrshed] = get_regmax(frame,blobsize,threshold)

sz = size(frame);
filt = conj(fftn(make_blob(sz,blobsize)));
smoothed_data = real(ifftn(fftn(frame).*filt));
regmax = smoothed_data.*(abs(smoothed_data) > threshold & imregionalmax(abs(smoothed_data)));
if nargout == 2
    wtrshed = watershed(-abs(smoothed_data)).*uint8(abs(smoothed_data)>threshold);
end

% remove regional maxima on the border
regmax(1,:,:)     = 0;
regmax(sz(1),:,:) = 0;
regmax(:,1,:)     = 0;
regmax(:,sz(2),:) = 0;