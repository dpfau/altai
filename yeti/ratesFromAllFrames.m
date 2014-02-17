load ROI_results.mat
data = padarray(loadframe(t),[0,0,1]);
rates = ratesFromFrame(data,ROIShapes,ROIOffset,numROI);
save(['rates/rates_' num2str(t)],'rates');