function ROI = extract_roi_morph(data,params)
% A morphological, rather than statistical, approach to identifying ROIs

sz = size(data);
ROI = [];
params.T = sz(end);
for i = 1:sz(end)
    fprintf('\nAnalyzing frame %d of %d',i,sz(end));
    ROI = update_roi_morph(data(:,:,i),ROI,params,i);
end
fprintf('\n')
% Post-processing steps

% Discard too-small ROIs
ROI = ROI(arrayfun(@(x)nnz(x.shape)>params.min_size,ROI));

% Separate accidentally overlapped ROIs with ICA
ROI_mat = cell2mat(arrayfun(@(x)x.shape(:).*x.prec(:),ROI,'UniformOutput',0));
[demixed_ROI,A,~] = fastica(ROI_mat');
demixed_ROI = diag(sign(mean(demixed_ROI,2)))*demixed_ROI; % Flip back any ROI that was flipped upside down
[~,idx] = max(abs(A')); % Find the index for the corresponding original ROI, since ICA permutes things
demixed_ROI = demixed_ROI(idx,:); % unpermute the ROIs
for i = 1:length(ROI)
    ROI(i).demixed = reshape(demixed_ROI(i,:),sz(1),sz(2)).*logical(ROI(i).shape);
    ROI(i).demixed = ROI(i).demixed/norm(ROI(i).demixed);
end

% % Split any ROIs with multiple maxima
% filt = conj(fft2(make_blob(sz(1),sz(2),params.blobsize)));
% retain = true(length(ROI),1);
% for i = 1:length(ROI)
%     smoothed_ROI = ifftn(fftn(ROI(i).demixed).*filt);
%     smoothed_ROI = smoothed_ROI .* (smoothed_ROI>0.1*max(smoothed_ROI(:)));
%     regmax = imregionalmax(smoothed_ROI);
%     if nnz(regmax) > 1 % We accidentally merged two blobs into one ROI
%         retain(i) = false;
%         wtrshed = watershed(-smoothed_ROI).*uint8(smoothed_ROI~=0);
%         idx = find(regmax);
%         for j = 1:length(idx)
%             region = wtrshed==wtrshed(idx(j));
%             if nnz(region) > params.min_size
%                 ROI(end+1).intensity = ROI(i).intensity;
%                 ROI(end).prec = ROI(i).prec .* region;
%                 ROI(end).shape = ROI(i).shape .* region;
%                 [x,y] = ind2sub(size(ROI(end).shape),idx(j));
%                 ROI(end).pos = [x,y];
%                 ROI(end).demixed = ROI(i).demixed .* region;
%             end
%         end
%     end
% end
% ROI = ROI(retain);