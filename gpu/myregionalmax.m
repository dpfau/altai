function regmax = myregionalmax(I)
% An easy regional maxima function. 
% Only looks at points where the value is greater than zero,
% and ignores the edges, but for our purposes it gives exactly 
% the same answer as Matlab's imregionalmax

switch ndims(I)
    case 2
        regmax = I(2:end-1,2:end-1)>0;
    case 3
        regmax = I(2:end-1,2:end-1,2:end-1)>0;
    otherwise
        error('Getting too fancy there, boyo. Stick to 2 or 3 dimensional arrays')
end

for i = -1:1
    for j = -1:1
        switch ndims(I)
            case 2
                if any([i,j])
                    regmax = regmax & ...
                        I(2:end-1,2:end-1) >= I(i+(2:end-1),j+(2:end-1));
                end
            case 3
                for k = -1:1
                    if any([i,j,k])
                        regmax = regmax & ...
                            I(2:end-1,2:end-1,2:end-1) >= I(i+(2:end-1),j+(2:end-1),k+(2:end-1));
                    end
                end
        end
    end
end
regmax = padarray(regmax,ones(1,ndims(I)));