function [fhExtractRatio] = ExtractRatio(vnChannels, bUsedRR)

% ExtractRatio - Construct an extraction function which uses the ratio of two channels
%
% Usage: [fhExtractRatio] = ExtractRatio(<vnChannels, bUsedRR>)
%
% 'vnChannels' is a two-element vector, which defines which channels to use
% in the ratio.  The respose will be channel 1 / channel 2.  The default is
% [1 2].
%
% 'bUsedRR' is a boolean flag (default false) which specifies whether to
% calculate delta R / R (if true).  If false, the raw ratio will be
% extracted.
%
% 'fhExtractRatio' will be a function handle that can be passed to
% ExtractRegionResponses, etc.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 28th October, 2011

% -- Default arguments

DEF_vnChannels = 1:2;
DEF_bUsedRR = false;


% -- Check arguments

if (~exist('vnChannels', 'var') || isempty(vnChannels))
   vnChannels = DEF_vnChannels;
end

if (~exist('bUsedRR', 'var') || isempty(bUsedRR))
   bUsedRR = DEF_bUsedRR;
end


% -- Return function handle

fhExtractRatio = @(fsData, vnPixels, vnFrames)fhExtractRatioFun(fsData, vnPixels, vnFrames, vnChannels, bUsedRR);


% --- END of ExtractRatio FUNCTION ---

   function [mfRawTrace, vfRegionTrace, fRegionResponse, nFramesInSample, vfPixelResponse] = fhExtractRatioFun(fsData, vnPixels, vnFrames, vnChannels, bUsedRR)
      
      tfTrace = double(fsData(vnPixels, vnFrames, vnChannels));
      
      % - Filter trace to remove zeros
      tfTrace(tfTrace == 0) = nan;
      
      mfRawTrace = tfTrace(:, :, 1) ./ tfTrace(:, :, 2);
      
      % - Calculate deltaR/R
      if (bUsedRR)
         vfBlankTrace = nanmean(double(fsData.BlankFrames(vnPixels, vnFrames)), 1);
         mfRawTrace = mfRawTrace ./ repmat(vfBlankTrace, size(mfRawTrace, 1), 1) - 1;
      end
      
      vfRegionTrace = nanmean(mfRawTrace, 1);
      fRegionResponse = nanmean(vfRegionTrace);
      vfPixelResponse = nanmean(mfRawTrace, 2);
      nFramesInSample = numel(vfRegionTrace);
   end
end
