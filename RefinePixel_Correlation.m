function sROIs = RefinePixel_Correlation(fsData, seedpixels, nChannel, fAlpha)
% finds pixel correlations to define ROIs
% based on RefinROIs_Correlation by D.Muir
% sROIs = ROI structure (see BWConComp for structure of sROIs
% fsData = output from FocusStack
% seedpixels = seed pixels in xy coordinates for which correlations are
%               calculated, xy is in columns, each row is a pixel
% nChannel = channel 1 for green, channel 2 for red
% fAlpha = ??
%
% 120621 by kampa@hifo.uzh.ch
%


% changed to read pixels instead of sROIs
for npixel = 1:size(seedpixels,1)
   % - Find the centre of mass of this ROI
   [vfI, vfJ] = ind2sub(size(fsData, 1:2), sROIs.PixelIdxList{npixel});
   fI = round(mean(vfI));
   fJ = round(mean(vfJ));
   
   % - Record the central pixel
   vnROIPixel(nROI) = sub2ind(size(fsData, 1:2), fI, fJ);
end

% - Get mean stack response
mfMeanAct = fsData.SummedAlignedFrames(:, :, :, nChannel) ./ size(fsData, 3);


% - Extract correlations
vnWindowStep = 1:20;
vnWindow = vnWindowStep(vnWindowStep <= sROIs.NumObjects);

tic;
tStart = toc;
fprintf(1, 'Measuring correlations: [%3d]%% %5d mins', 0, 0);
while (~isempty(vnWindow))
   % - Extract activity traces
   mfROIActivity = fsData.AlignedStack(vnROIPixel(vnWindow), :, nChannel);
   
   % - Measure correlations
   tfCorr(:, :, vnWindow) = CorrelationImage(fsData, mfROIActivity, 1:size(fsData, 3), nChannel, mfMeanAct);
   
   % - Display some progress
   tCurr = toc;
   fPropDone = vnWindow(end) / sROIs.NumObjects;
   fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3d]%% %5d mins', ...
      floor(fPropDone*100), ceil((tCurr - tStart) / fPropDone * (1-fPropDone)/60));
  
   % - Move to next window
   vnWindow = vnWindow(end) + vnWindowStep;
   vnWindow = vnWindow(vnWindow <= sROIs.NumObjects);
end

% % - Find the correlation between this point and the rest of the imaged region
% tfCorr(:, :, nROI) = CorrelationImage(fsData, {fI, fJ, ':', nChannel});
% 

% - Find ROIs using this threshold
for (nROI = 1:sROIs.NumObjects)
   % - Determine a threshold for this correlation
   vfCorrs = sort(reshape(tfCorr(:, :, nROI), 1, []));
   fCorrThreshold = vfCorrs(ceil(numel(vfCorrs) * (1-fAlpha)));
   
   % - Find pixels with high correlation
   mbInclude = tfCorr(:, :, nROI) > fCorrThreshold;
   
   % - Find the centre of mass of this ROI
   [vfI, vfJ] = ind2sub(size(fsData, 1:2), sROIs.PixelIdxList{nROI});
   fI = round(mean(vfI));
   fJ = round(mean(vfJ));

   % - Include only regions connected to centre of mass pixel
   sThisRoi = bwconncomp(bwselect(mbInclude, fJ, fI, 4));
   sROIs.PixelIdxList{nROI} = sThisRoi.PixelIdxList{1};
end


% --- END of RefineROIs_Correlation.m ---
