function make_regmax_vid(data,regmax,z)

figure
sz = size(data);
assert(z>1&z<sz(3));
img = mat2gray(data(:,:,z-1:z+1,:));
foo = [' '*ones(1,max(0,z-3)) 'rrgbb'];
if z == 2, foo = foo(2:end); end
vidObj = VideoWriter(['RegionalMaxima' num2str(z) '.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 24;
open(vidObj);
for i = 1:sz(4)
    image(img(:,:,:,i));
    axis image
    hold on
    for j = 5:9
        w = find(regmax(:,:,j,i));
        [x,y] = ind2sub(sz(1:2),w);
        k = length(x);
        scatter(y,x,2*abs(regmax(sub2ind(sz,x,y,j*ones(k,1),i*ones(k,1)))),[foo(j) 'o'])
    end
    hold off
    drawnow
    writeVideo(vidObj, getframe(gca));
end
close(vidObj)