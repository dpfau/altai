function ROI = update_roi_morph(frame,ROI,params,t)

[regmax,wtrshed] = get_regmax(frame,...
                    params.blobsize,...
                    params.threshold);
          
[x,y] = ind2sub(size(regmax),find(regmax));
nMax = length(x); % number of regional maxima
wtrshed_id = zeros(nMax,1);
for i = 1:nMax % match watershed labels with regional maxima
    wtrshed_id(i) = wtrshed(x(i),y(i));
end
nROI = length(ROI); % number of ROI from previous steps
fprintf('...Found %d regional maxima',nMax);

overlap = false(nMax,nROI);
distance = zeros(nMax,nROI);
for i = 1:nMax % Find which existing ROI fall within a watershed
    region = wtrshed==wtrshed_id(i);
    for j = 1:nROI
        distance(i,j) = norm(ROI(j).pos-[x(i),y(i)]);
        overlap(i,j) = any(region(:)&ROI(j).shape(:));
    end
end
new_roi = zeros(nMax,1);
assignment = nROI+(1:nMax);
intensity = regmax(logical(regmax)); % Height of the regional maxima.
if nROI > 0
    ROI_mat = cell2mat(arrayfun(@(x)x.shape(:),ROI,'UniformOutput',0));
    rates = pinv(ROI_mat)*frame(:);
    resid_var = params.var_offset*ones(numel(frame),1);
    for j = 1:nROI
        resid_var(logical(ROI(j).prec)) = resid_var(logical(ROI(j).prec)) + params.var_slope*rates(j) + rates(j)^2./ROI(j).prec(logical(ROI(j).prec));
    end
    resid = (frame(:) - ROI_mat*rates)./sqrt(resid_var);
else
    resid = frame;
end
for i = 1:nMax
    j = 0; % Index of the ROI to be merged into, or 0 if there is nothing to merge into
%     if any(distance(i,:) < params.regmax_dist)
%         [~,j] = min(distance(i,:));
%         assignment(i) = j;
    if any(overlap(i,:))
        region = wtrshed(:)==wtrshed_id(i)&(sum(ROI_mat(:,overlap(i,:)),2)); % The region inside the watershed that is already included in other ROIs should not be significantly different from background noise, unless there is a new ROI in this watershed
        if 1 - chi2cdf(norm(resid(region))^2,nnz(region)) > params.pval % If it is not significantly different, we have to pick an ROI to merge it with. Let's go with the one that has the greatest power inside the region of overlap.
            idx = find(overlap(i,:));
            pow = zeros(length(idx),1);
            for ii = 1:length(idx)
                pow(ii) = norm(rates(idx(ii))*ROI_mat(region,idx(ii)));
            end
            [~,jj] = max(pow);
            j = idx(jj);
        end
    end
    % This whole section still needs to be updated to reflect the possible
    % influence of multiple ROIs in one watershed
    if j ~= 0
        ROI(j).shape = (ROI(j).prec .* ROI(j).shape + ...
            intensity(i)*(wtrshed==wtrshed_id(i)) .* resid /(params.var_offset + params.var_slope*intensity(i))) ./ ...
            (ROI(j).prec + intensity(i)^2*(wtrshed==wtrshed_id(i))/(params.var_offset + params.var_slope*intensity(i)));
        ROI(j).shape(isnan(ROI(j).shape)) = 0;
        ROI(j).pos = (sum(ROI(j).intensity.^2) * ROI(j).pos + intensity(i).^2 * [x(i),y(i)])...
            /(sum(ROI(j).intensity.^2) + intensity(i).^2);
        ROI(j).intensity(t) = intensity(i);
        ROI(j).prec = ROI(j).prec + intensity(i)^2*(wtrshed==wtrshed_id(i))/(params.var_offset + params.var_slope*intensity(i));
    else
        new_roi(i) = 1;
        region = wtrshed==wtrshed_id(i);
        ROI(end+1).intensity = zeros(params.T,1);
        ROI(end).intensity(t) = intensity(i);
        ROI(end).prec = intensity(i)^2*(wtrshed==wtrshed_id(i))/(params.var_offset + params.var_slope*intensity(i));
        ROI(end).shape = frame.*region/intensity(i);
        ROI(end).pos = [x(i),y(i)];
    end
end

fprintf('...Added %d new regions of interest\n',nnz(new_roi));

if isfield(params,'optose') && params.optose == true
    h = zeros(nMax,1); % handles to the text objects in each frame
    figure(61)
    colormap gray
    hold off
    imshow(mat2gray(frame))
    axis image
    hold on
    himg = imagesc(label2rgb(wtrshed));
    set(himg,'AlphaData',.25);
    scatter(y,x,50,'ro','LineWidth',2)
    
    for i = 1:length(ROI)
        scatter(ROI(i).pos(2),ROI(i).pos(1),100,'bo','LineWidth',2);
    end
    drawnow
end