function [fhExtractMean] = ExtractMean(nChannel, bUsedFF)

% ExtractMean - FUNCTION Extract the time and space average of a response trace
%
% Usage: [fhExtractMean] = ExtractMean(<nChannel, bUsedFF>)
%
% This function extracts the time and space average of an ROI response
% trace, optionally computing the delta F / F values.
%
% 'nChannel' specifies which channel to extract data from. 'bUsedFF' is a
% boolean flag specifying whether or not to extract delta F / F values
% using the pre-defined blanks (default: off).

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 3rd November, 2011

% -- Default arguments

DEF_nChannel = 1;
DEF_bUsedFF = false;


% -- Check arguments

if (~exist('nChannel', 'var') || isempty(nChannel))
   nChannel = DEF_nChannel;
end

if (~exist('bUsedFF', 'var') || isempty(bUsedFF))
   bUsedFF = DEF_bUsedFF;
end

% -- Return function handle

fhExtractMean = @(fsData, vnPixels, vnFrames)fhExtractMeanFun(fsData, vnPixels, vnFrames, nChannel, bUsedFF);

% --- END of ExtractMean FUNCTION ---

   function [mfRawTrace, vfRegionTrace, fRegionResponse, nFramesInSample, vfPixelResponse] = ...
         fhExtractMeanFun(fsData, vnPixels, vnFrames, nChannel, bUsedFF)
      
      mfRawTrace = double(fsData(vnPixels, vnFrames, nChannel));      
      
      % - Calculate deltaF/F
      if (bUsedFF)
         vfBlankTrace = nanmean(double(fsData.BlankFrames(vnPixels, vnFrames)), 1);
         mfRawTrace = mfRawTrace ./ repmat(vfBlankTrace, size(mfRawTrace, 1), 1) - 1;
      end
      
      vfRegionTrace = nanmean(mfRawTrace, 1);
      fRegionResponse = nanmean(vfRegionTrace);
      vfPixelResponse = nanmean(mfRawTrace, 2);
      nFramesInSample = numel(vfRegionTrace);
   end
end

% --- END of ExtractMean.m ---
