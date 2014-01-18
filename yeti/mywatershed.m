function idx = mywatershed(I,seed)
% An algorithn to compute one watershed basin starting from a seed. When
% used to create a closure and passed to arrayfun with a gpuArray as the
% second argument, is automatically executed in parallel on the GPU.

n = ndims(I);
sz = size(I);

idx = seed;
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
    for i = -1:1
        for j = -1:1
            switch n
                case 2
                    if any([i,j])
                        idx_shift = setdiff(sub2ind(sz,min(sz(1),max(1,x+i)),min(sz(2),max(1,y+j))),idx);
                        boundary = grow_boundary(boundary, idx_shift, I);
                    end
                case 3
                    for k = -1:1
                        if any([i,j,k])
                            idx_shift = setdiff(sub2ind(sz,min(sz(1),max(1,x+i)),min(sz(2),max(1,y+j)),min(sz(3),max(1,z+k))),idx);
                            boundary = grow_boundary(boundary, idx_shift, I);
                        end
                    end
            end
        end
    end
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
    idx = [idx; boundary(in_watershed)];
    x = [x; xb(in_watershed)];
    y = [y; yb(in_watershed)];
    if n == 3
        z = [z; zb(in_watershed)];
    end
    npix = nnz(in_watershed);
end

function boundary = grow_boundary(boundary, idx_shift, I)

nonzero   = I(idx_shift)>0;
idx_shift = idx_shift(nonzero);
overlap   = ismember(idx_shift,boundary);
boundary  = [boundary; idx_shift(~overlap)];