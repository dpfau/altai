load ROI_results.mat
data = loadframe(t);
rates = ratesFromFrame(data,ROIShapes,ROIOffet,numROI);
save(['rates_' num2str(t)],'rates');