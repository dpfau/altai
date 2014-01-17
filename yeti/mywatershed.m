function watershed = mywatershed(I,marker)

idx = find(marker); % index of the pixels in all watersheds
label = (1:length(idx))'; % watershed label for the particular pixel
n = ndims(I);
sz = size(I);
switch n
    case 2
        [x,y] = ind2sub(sz,idx);
    case 3
        [x,y,z] = ind2sub(sz,idx);
    otherwise
        error('Hold your horses, I only implemented this for 2D and 3D')
end

npix = 1; % number of pixels added to watershed in the last iteration
while npix > 0
    % compute boundary of watershed
    boundary       = [];
    boundary_label = [];
    overlap_different = [];
    for i = -1:1
        for j = -1:1
            switch n
                case 2
                    if any([i,j])
                        [idx_shift, ia] = setdiff(sub2ind(sz,min(sz(1),max(1,x+i)),min(sz(2),max(1,y+j))),idx);
                        [boundary, boundary_label, overlap_different] = grow_boundary(boundary, boundary_label, overlap_different, idx_shift, ia, I, label);
                    end
                case 3
                    for k = -1:1
                        if any([i,j,k])
                            [idx_shift, ia] = setdiff(sub2ind(sz,min(sz(1),max(1,x+i)),min(sz(2),max(1,y+j)),min(sz(3),max(1,z+k))),idx);
                            [boundary, boundary_label, overlap_different] = grow_boundary(boundary, boundary_label, overlap_different, idx_shift, ia, I, label);
                        end
                    end
            end
        end
    end
    [boundary,ia] = setdiff(boundary, overlap_different);
    boundary_label = boundary_label(ia);
    switch n
        case 2
            [xb,yb]    = ind2sub(sz,boundary);
        case 3
            [xb,yb,zb] = ind2sub(sz,boundary);
    end
    
    % find the index and value of the maximum neighbor of each boundary pixel
    neighbor_val = zeros(size(boundary));
    neighbor_idx = zeros(size(boundary));
    for i = -1:1
        for j = -1:1
            switch ndims(I)
                case 2
                    if any([i,j])
                        subs = sub2ind(sz,min(sz(1),max(1,xb+i)),min(sz(2),max(1,yb+j)));
                        greater = I(subs) > neighbor_val;
                        neighbor_val(greater) = I(subs(greater));
                        neighbor_idx(greater) = subs(greater);
                    end
                case 3
                    for k = -1:1
                        if any([i,j,k])
                            subs = sub2ind(sz,min(sz(1),max(1,xb+i)),min(sz(2),max(1,yb+j)),min(sz(3),max(1,zb+k)));
                            greater = I(subs) > neighbor_val;
                            neighbor_val(greater) = I(subs(greater));
                            neighbor_idx(greater) = subs(greater);
                        end
                    end
                otherwise
                    error('Whoa there, I only implemented this for 2D and 3D');
            end
        end
    end
    
    % if the maximum value is inside this watershed already, add the pixel to the watershed
    [in_watershed] = ismember(neighbor_idx,idx);
    idx   = [idx;   boundary(in_watershed)];
    label = [label; boundary_label(in_watershed)];
    x = [x; xb(in_watershed)];
    y = [y; yb(in_watershed)];
    if n == 3
        z = [z; zb(in_watershed)];
    end
    npix = nnz(in_watershed);
    fprintf('.')
end
fprintf('\n')
watershed = sparse(idx,ones(size(idx)),label,numel(I),1);

function [boundary, boundary_label, overlap_different] = grow_boundary(boundary, boundary_label, overlap_different, idx_shift, ia, I, label)

label_shift = label(ia);

nonzero = I(idx_shift)>0;
idx_shift   = idx_shift(nonzero);
label_shift = label_shift(nonzero);

[overlap, loc] = ismember(idx_shift,boundary);
overlap_different_idx = overlap;
overlap_different_idx(overlap) = boundary_label(loc(overlap)) ~= label_shift(overlap);
overlap_different = [overlap_different; idx_shift(overlap_different_idx)];

boundary       = [boundary;       idx_shift(~overlap)];
boundary_label = [boundary_label; label_shift(~overlap)];