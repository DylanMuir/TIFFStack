function [fhExtractionFun] = ExtractPeak(nSampleFrames, nChannel, bUsedFF)

% ExtractPeak - Extraction function for returning the peak of the response
%
% Usage: [fhExtractionFun] = ExtractPeak(<nSampleFrames, nChannel, bUsedFF>)
%
% This function will extract the trace of a region, find the peak point
% during that trace, and return the average of a few frames around that
% point.
%
% 'nSampleFrames' specifies the number of frames to sample for the mean of
% the peak (default: 5);
%
% 'nChannel' specifies which channel to extract (default: 1).
%
% 'bUsedFF' specifies whether or not to extract the delta F / F value (if
% true), or the raw signal response (if false)
%
% 'fhExtractionFun' will be a function handle that can be passed to
% ExtractRegionResponses, etc.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 29th October, 2011

% -- Defaults

DEF_nSampleFrames = 5;
DEF_nChannel = 1;
DEF_bUsedFF = false;


% -- Check arguments

if (~exist('nSampleFrames', 'var') || isempty(nSampleFrames))
   nSampleFrames = DEF_nSampleFrames;
end

if (~exist('nChannel', 'var') || isempty(nChannel))
   nChannel = DEF_nChannel;
end

if (~exist('bUsedFF', 'var') || isempty(bUsedFF))
   bUsedFF = DEF_bUsedFF;
end


% -- Return function handle

fhExtractionFun = @(fsData, vnPixels, vnFrames)fhExtractPeak(fsData, vnPixels, vnFrames, nSampleFrames, nChannel, bUsedFF);


% --- END of ExtractPeak.m ---

   function [mfRawTrace, vfRegionTrace, fRegionResponse, nFramesInSample, vfPixelReponse] = fhExtractPeak(fsData, vnPixels, vnFrames, nSampleFrames, nChannel, bUsedFF)
      
      mfRawTrace = double(fsData(vnPixels, vnFrames, nChannel));
      
      % - Calculate deltaF/F
      if (bUsedFF)
         vfRegionBlankTrace = nanmean(fsData.BlankFrames(vnPixels, vnFrames), 1);
         mfRawTrace = mfRawTrace ./ repmat(vfRegionBlankTrace, size(mfRawTrace, 1), 1) - 1;
      end
      
      vfRegionTrace = nanmean(mfRawTrace, 1);
      
      % - Locate peak
      [nul, nPeakIndex] = nanmax(vfRegionTrace);
      
      % - Average five frames surrounding peak
      vnWindow = (nPeakIndex - floor(nSampleFrames/2)):(nPeakIndex + ceil(nSampleFrames/2)-1);
      vnWindow = vnWindow(vnWindow >= 1);
      vnWindow = vnWindow(vnWindow <= numel(vfRegionTrace));
      
      fRegionResponse = nanmean(vfRegionTrace(vnWindow));
      nFramesInSample = numel(vnWindow);
      
      if (nargout > 4)
         vfPixelReponse = nanmax(mfRawTrace, [], 2);
      end
   end

end