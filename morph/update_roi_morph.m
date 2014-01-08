function ROI = update_roi_morph(frame,ROI,params,t)

[regmax,wtrshed] = get_regmax(frame,...
                    params.blobsize,...
                    params.threshold);
          
[x,y] = ind2sub(size(regmax),find(regmax));
k = length(x); % number of regional maxima
wtrshed_id = zeros(k,1);
for i = 1:k % match watershed labels with regional maxima
    wtrshed_id(i) = wtrshed(x(i),y(i));
end
q = length(ROI); % number of ROI from previous steps
fprintf('...Found %d regional maxima',k);

overlap = zeros(k,q,2);
distance = zeros(k,q);
for i = 1:k % Find which existing ROI fall within a watershed
    pix = wtrshed==wtrshed_id(i);
    for j = 1:q
        distance(i,j) = norm(ROI(j).pos-[x(i),y(i)]);
        overlap(i,j,1) = nnz(pix&ROI(j).shape)/nnz(pix);
        overlap(i,j,2) = nnz(pix&ROI(j).shape)/nnz(ROI(j).shape);
    end
end
new_roi = zeros(k,1);
assignment = q+(1:k);
intensity = regmax(logical(regmax));
for i = 1:k
    j = 0;
    if any(distance(i,:) < params.regmax_dist) || any(overlap(i,:,1) > 0.5) || any(overlap(i,:,2) > 0.5)
        j = find(distance(i,:) < 3);
        if length(j) > 1
            [~,j] = min(distance(i,:));
        elseif isempty(j)
            idx = find(overlap(i,:,1)>0.5 | overlap(i,:,2)>0.5);
            pval = zeros(length(idx),1);
            for j = 1:length(idx)
                region = ROI(idx(j)).shape & wtrshed==wtrshed_id(i);
                resid = (frame(region)/intensity(i) - ROI(idx(j)).shape(region))./sqrt(1./(params.prec*intensity(i)^2) + 1./ROI(idx(j)).prec(region)); % residual scaled so each pixel should be i.i.d. unit normal
                pval(j) = 1 - chi2cdf(norm(resid)^2,nnz(region));
            end
            j = find(pval>params.pval);
            if length(j) > 1 || isempty(j) % If you can't reject anything, it's probably so ambiguous that it's better to toss it out
                j = 0;
            else
                j = idx(j);
            end
        end
        assignment(i) = j;
    end
    if j ~= 0
        ROI(j).shape = (ROI(j).prec .* ROI(j).shape + ...
            params.prec*intensity(i)*(wtrshed==wtrshed_id(i)) .* frame) ./ ...
            (ROI(j).prec + params.prec*intensity(i)^2*(wtrshed==wtrshed_id(i)));
        ROI(j).shape(isnan(ROI(j).shape)) = 0;
        ROI(j).pos = (sum(ROI(j).intensity.^2) * ROI(j).pos + intensity(i).^2 * [x(i),y(i)])...
            /(sum(ROI(j).intensity.^2) + intensity(i).^2);
        ROI(j).intensity(t) = intensity(i);
        ROI(j).prec = ROI(j).prec + params.prec*intensity(i)^2*(wtrshed==wtrshed_id(i));
    else
        new_roi(i) = 1;
        ROI(end+1).intensity = zeros(params.T,1);
        ROI(end).intensity(t) = intensity(i);
        ROI(end).prec = params.prec*intensity(i)^2*(wtrshed==wtrshed_id(i));
        ROI(end).shape = frame.*(wtrshed==wtrshed_id(i))/intensity(i);
        ROI(end).pos = [x(i),y(i)];
    end
end

% connected = sparse(length(ROI),length(ROI));
% for i = 1:length(ROI)
%     for j = i+1:length(ROI)
%         connected(i,j) = norm(ROI(i).pos-ROI(j).pos)<3 ;%|| ...
%             %(ROI(i).shape(floor(ROI(j).pos(1)),floor(ROI(j).pos(2))) ~= 0 && ...
%             % ROI(j).shape(floor(ROI(i).pos(1)),floor(ROI(i).pos(2))) ~= 0);
%     end
% end
% connected = connected + connected';
% connected(1:length(ROI)+1:end) = 1;
% [S,C] = graphconncomp(connected);
% for i = 1:S
%     ROI_(i).intensity = zeros(params.T,1);
%     ROI_(i).shape = zeros(size(frame));
%     ROI_(i).pos = [0,0];
%     for j = 1:length(ROI)
%         if C(j)==i
%             ROI_(i).intensity = ROI_(i).intensity + ROI(j).intensity;
%             ROI_(i).shape = ROI_(i).shape + sum(ROI(j).intensity) * ROI(j).shape;
%             ROI_(i).pos = ROI_(i).pos + sum(ROI(j).intensity) * ROI(j).pos;
%         end
%     end
%     ROI_(i).shape = ROI_(i).shape / sum(ROI_(i).intensity);
%     ROI_(i).pos = ROI_(i).pos / sum(ROI_(i).intensity);
% end
% 
% fprintf('...Added %d new regions of interest, merged %d regions of interest\n',nnz(new_roi),length(ROI)-length(ROI_));
% ROI = ROI_;
fprintf('...Added %d new regions of interest\n',nnz(new_roi));

if isfield(params,'optose') && params.optose == true
    h = zeros(k,1); % handles to the text objects in each frame
    figure(61)
    colormap gray
    hold off
    imshow(mat2gray(frame))
    axis image
    hold on
    himg = imagesc(label2rgb(wtrshed));
    set(himg,'AlphaData',.25);
    scatter(y,x,50,'ro','LineWidth',2)
%     for i = 1:k
%         if assignment(i) ~= 0
%             h(i) = text(y(i)+10,x(i),num2str(assignment(i)),'Color','k','FontSize',14,'LineWidth',1.5);
%         end
%     end
    
    for i = 1:length(ROI)
        % roi_edge = bwboundaries(ROI(i).shape~=0);
        % line([roi_edge{1}(:,2),[roi_edge{1}(2:end,2);roi_edge{1}(1,2)]],[roi_edge{1}(:,1),[roi_edge{1}(2:end,1);roi_edge{1}(1,1)]],'Color','r');
        scatter(ROI(i).pos(2),ROI(i).pos(1),100,'bo','LineWidth',2);
    end
    drawnow
end