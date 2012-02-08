function PlotSpatialResponse(sAnalysis, vnStimuliSum)

% - Make a blank image
mfImage = zeros(sAnalysis.sCellRegions.ImageSize);

for (nRegion = 1:sAnalysis.sCellRegions.NumObjects)
   mfImage(sAnalysis.sCellRegions.PixelIdxList{nRegion}) = sum(sAnalysis.mfStimMeans(nRegion, vnStimuliSum));
end

figure;
imagesc(mfImage);
colorbar;
cMap = jet(256);
vCAxis = caxis;
caxis([0 vCAxis(2)]);

cMap(1, :) = 0;
colormap(cMap);
