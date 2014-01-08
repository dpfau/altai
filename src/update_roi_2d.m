function ROI = update_roi_2d(frame,ROI,params)

sz = size(frame);
[regmax,wtrshed] = get_regmax(frame,...
                    params.blobsize,...
                    params.threshold);
          
[x,y] = ind2sub(sz,find(regmax));
k = length(x); % number of regional maxima
q = length(ROI); % number of ROI from previous steps
fprintf('...Found %d regional maxima',k);
if isfield(params,'optose') && params.optose == true
    h = zeros(k,1); % handles to the text objects in each frame
    figure(61)
    subplot(1,2,1)
    colormap gray
    hold off
    imshow(mat2gray(frame))
    axis image
    hold on
    himg = imagesc(label2rgb(wtrshed));
    set(himg,'AlphaData',.15);
    scatter(y,x,50,'g.')
    for i = 1:k
        h(i) = text(y(i)+10,x(i),num2str(q+i),'Color','g');
    end
    
    for i = 1:q
        scatter(ROI(i).pos(2),ROI(i).pos(1),100,'bo');
    end
end

% First, add all regmax in as potential ROI
[mu,Sig] = make_roi_prior(params.order);
for i = 1:k
    ROI(q+i).pos = [x(i),y(i)];
    ROI(q+i).mu = mu;
    ROI(q+i).Sig = Sig;
    ROI(q+i).nFrames = 0; % number of frames the ROI appears in
    ROI(q+i).maxIntensity = 0;
end

distance = zeros(k,q+k); % technically, the \ell_\infty norm between two ROI locations. If greater than the window size, there is no overlap.
% SLOW! Replace with KD-tree to scale?
for i = 1:k % Find which regional maxima are close enough to existing ROI that they might be the same
    for j = 1:q+k
        distance(i,j) = max(abs([x(i),y(i)]-ROI(j).pos));
    end
end
assignment = q+(1:k); % Initially, each regmax is assigned to the corresponding new ROI
inverse_assignment = [zeros(1,q),1:k]; % Really the style that MATLAB forces you to write in is so kludgey compared to having pointers.

% Try the trivial thing of merging any regmaxes that are close enough
% together. I know this won't work on all data, but I want to see how
% common failures are.
neighbors = distance <= 5;
for i = 1:k
    if any(neighbors(i,1:q))
        assignment(i) = find(neighbors(i,1:q),1);
        inverse_assignment(q+i) = 0;
        inverse_assignment(find(neighbors(i,1:q),1)) = i;
        if isfield(params,'optose') && params.optose == true
            set(h(i),'String',num2str(assignment(i)))
        end
    end
end

% Then, Gibbs sample regional maxima, figuring out which ROI they can be assigned to.
intensity = regmax(logical(regmax));
% delta = ones(k,1);
% while any(delta)
%     delta = zeros(k,1);
%     for i = 1:k
%         idx = find(distance(i,:) < params.window/2 & inverse_assignment == 0); % the threshold for how far a regmax can be from an ROI is a little stricter than the window size
%         if ~isempty(idx)
%             j_ = assignment(i);
%             neighbors = distance(i,:) < params.window & inverse_assignment ~=0; % We only care about ROIs that are assigned to a regmax. We assume for simplicity that all other ROIs are silent. Sometimes this is wrong but hopefully it won't trip us up too much.
%             loglik_ = get_evidence(frame,ROI(neighbors),params.basis,intensity(inverse_assignment(neighbors)),params.sig,[x(i),y(i)]);
%             logliks = zeros(length(idx),1);
%             inverse_assignment(j_) = 0;
%             neighbors(j_) = 0;
%             for j = 1:length(idx)
%                 inverse_assignment(idx(j)) = i;
%                 neighbors(idx(j)) = 1;
%                 logliks(j) = get_evidence(frame,ROI(neighbors),params.basis,intensity(inverse_assignment(neighbors)),params.sig,[x(i),y(i)]);
%                 inverse_assignment(idx(j)) = 0;
%                 neighbors(idx(j)) = 0;
%             end
%             j = sample_logliks([logliks;loglik_]);
%             if j ~= length(idx)+1
%                 delta(i) = 1;
%                 assignment(i) = idx(j);
%                 inverse_assignment(idx(j)) = i;
%                 if isfield(params,'optose') && params.optose == true
%                     set(h(i),'String',num2str(idx(j)))
%                     drawnow
%                 end
%             else
%                 inverse_assignment(j_) = i;
%             end
%         end
%     end
% end
fprintf('...Added %d new regions of interest.',nnz(inverse_assignment(q+1:end)));

% Then, having identified which regmax's are new ROIs and which are existing ones, update the shape given the new data
ROI(assignment) = update_posterior(frame,ROI(assignment),params.basis,intensity,params.sig);

if isfield(params,'optose') && params.optose == true
    resid = frame; % residual
    foo = zeros(size(frame));
    [m,n,p] = size(params.basis);
    for i = 1:k
        roi_shape = reshape(reshape(params.basis,m*n,p)*ROI(assignment(i)).mu,m,n);
        offset = ROI(assignment(i)).pos - floor([m,n]/2);
        foo(offset(1)+(1:m),offset(2)+(1:n)) = roi_shape;
        resid = resid - intensity(i)*foo;
        foo(offset(1)+(1:m),offset(2)+(1:n)) = 0;
    end
    figure(61)
    subplot(1,2,2)
    colormap gray
    imshow(mat2gray(resid))
    axis image
end

% And remove new ROIs that aren't necessary
ROI = ROI([true(1,q),inverse_assignment(q+1:end)~=0]);

if isfield(params,'optose') && params.optose == true
    figure(61)
    subplot(1,2,1)
    for i = q+1:length(ROI)
        scatter(ROI(i).pos(2),ROI(i).pos(1),100,'ro');
    end
    drawnow
end

if isfield(params,'pause')  && params.pause == true, pause; end