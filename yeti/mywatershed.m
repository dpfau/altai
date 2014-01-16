function watershed = mywatershed(I,marker)

watershed = zeros(size(I),'uint16');
idx = find(marker);
for id = 1:length(idx)
    watershed(idx(id)) = id;
end

npix = 1; % number of pixels added to watershed in the last iteration
while npix > 0
    % compute boundary of watershed
    boundary = zeros(size(watershed),'uint16');
    overlap = false(size(watershed));
    overlap_different = false(size(watershed));
    for i = -1:1
        for j = -1:1
            switch ndims(I)
                case 2
                    if any([i,j])
                        delta = ( watershed(i+(2:end-1),j+(2:end-1)) .* uint16( watershed(2:end-1,2:end-1) == 0 ) );
                        overlap(2:end-1,2:end-1) = overlap(2:end-1,2:end-1) | ...
                            ( delta & boundary(2:end-1,2:end-1) );
                        overlap_different(2:end-1,2:end-1) = overlap_different(2:end-1,2:end-1) | ...
                            ( delta & boundary(2:end-1,2:end-1) & (boundary(2:end-1,2:end-1) ~= delta) );
                        boundary(2:end-1,2:end-1) = boundary(2:end-1,2:end-1) + ( delta .* uint16( ~overlap(2:end-1,2:end-1) ) );
                    end
                case 3
                    for k = -1:1
                        if any([i,j,k])
                            delta = ( watershed(i+(2:end-1),j+(2:end-1),k+(2:end-1)) .* uint16( watershed(2:end-1,2:end-1,2:end-1) == 0 ) );
                            overlap(2:end-1,2:end-1,2:end-1) = overlap(2:end-1,2:end-1,2:end-1) | ...
                                ( delta & boundary(2:end-1,2:end-1,2:end-1) );
                            overlap_different(2:end-1,2:end-1,2:end-1) = overlap_different(2:end-1,2:end-1,2:end-1) | ...
                                ( delta & boundary(2:end-1,2:end-1,2:end-1) & (boundary(2:end-1,2:end-1,2:end-1) ~= delta) );
                            boundary(2:end-1,2:end-1,2:end-1) = boundary(2:end-1,2:end-1,2:end-1) + ( delta .* uint16( ~overlap(2:end-1,2:end-1,2:end-1) ) );
                        end
                    end
                otherwise
                    error('Whoa there, I only implemented this for 2D and 3D');
            end
        end
    end
    boundary = boundary .* uint16( (I>0) & ~overlap_different );
    boundary_idx = find(boundary);
    switch ndims(I)
        case 2
            [x,y] = ind2sub(size(I),boundary_idx);
        case 3
            [x,y,z] = ind2sub(size(I),boundary_idx);
        otherwise
            error('Whoa there, I only implemented this for 2D and 3D');
    end
    
    % find the index and value of the maximum neighbor of each boundary pixel
    neighbor_val = zeros(size(boundary_idx));
    neighbor_idx = zeros(size(boundary_idx));
    for i = -1:1
        for j = -1:1
            switch ndims(I)
                case 2
                    if any([i,j])
                        subs = sub2ind(size(I),x+i,y+j);
                        greater = I(subs) > neighbor_val;
                        neighbor_val(greater) = I(subs(greater));
                        neighbor_idx(greater) = subs(greater);
                    end
                case 3
                    for k = -1:1
                        if any([i,j,k])
                            subs = sub2ind(size(I),x+i,y+j,z+k);
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
    in_watershed = watershed(neighbor_idx);
    watershed(boundary_idx) = in_watershed;
    npix = nnz(in_watershed);
    fprintf('.')
end
