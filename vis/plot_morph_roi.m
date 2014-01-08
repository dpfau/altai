function plot_morph_roi(ROI)

sz = size(ROI(1).shape);
img = zeros([sz,3]);
for i = 1:length(ROI)
    roi_img = cat(3,rand*ones(sz),mat2gray(ROI(i).shape),mat2gray(ROI(i).prec));
    img = img + hsv2rgb(roi_img);
end

figure(62)
clf
image(mat2gray(img));
axis image
hold on
for i = 1:length(ROI)
    scatter(ROI(i).pos(2),ROI(i).pos(1),100,'ro','LineWidth',2)
end