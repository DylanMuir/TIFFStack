function fhExtFun = dFFFun(nChannel)

   fhExtFun = @(fsData, vnPixels, vnFrames)fhExtractdFF(fsData, vnPixels, vnFrames, nChannel);

   function [mfRawTrace, vfRegionTrace, fRegionResponse, nFramesInSample, vfPixelResponse] = fhExtractdFF(fsData, vnPixels, vnFrames, nChannel)
      
      mfRawTrace = fsData(vnPixels, vnFrames, nChannel);
      mfBlankTrace = fsData.BlankFrames(vnPixels, vnFrames);
      mfRawTrace = mfRawTrace ./ mfBlankTrace - 1;
      vfRegionTrace = nanmean(mfRawTrace, 1);
      fRegionResponse = nanmean(vfRegionTrace);
      vfPixelResponse = nanmean(mfRawTrace, 2);
      nFramesInSample = numel(vfRegionTrace);
   end

end
