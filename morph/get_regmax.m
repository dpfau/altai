function [regmax,wtrshed] = get_regmax(frame,blobsize,threshold)

sz = size(frame);
if length(sz) == 2
    filt = conj(fft2(make_blob(sz(1),sz(2),blobsize)));
else
    filt1 = conj(fft2(make_blob(sz(1),sz(2),blobsize)));
    filt = zeros(sz(1:3));
    filt(:,:,1) = filt1;
    filt(:,:,2) = 0.25*filt1;
    filt(:,:,end) = 0.25*filt1;
    filt = fft(filt,[],3);
end

smoothed_data = ifftn(fftn(frame).*filt);
regmax = smoothed_data.*(abs(smoothed_data) > threshold & imregionalmax(abs(smoothed_data)));
if nargout == 2
    wtrshed = watershed(-abs(smoothed_data)).*uint8(abs(smoothed_data)>threshold);
end

% remove regional maxima on the border
regmax(1,:,:)     = 0;
regmax(sz(1),:,:) = 0;
regmax(:,1,:)     = 0;
regmax(:,sz(2),:) = 0;