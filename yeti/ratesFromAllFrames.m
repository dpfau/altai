load ROI_results.mat
data = loadframe(t);
rates = ratesFromFrame(data,ROIShapes,ROIOffset,numROI);
save(['rates_' num2str(t)],'rates');