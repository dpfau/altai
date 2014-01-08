function plot_regmax_frame(data,sub,z,t)

colormap gray
imagesc(data(:,:,z,t),[min(data(:)),max(data(:))])
hold on
for k = find(sub(:,4)==t&sub(:,3)==z)'
    scatter(sub(k,2),sub(k,1),50,'g.'); 
end
drawnow
hold off