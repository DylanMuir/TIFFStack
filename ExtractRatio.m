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

fhExtractRatio = @(fsData, cvnPixels, vnFrames)fhExtractRatioFun(fsData, cvnPixels, vnFrames, vnChannels, bUsedRR);


% --- END of ExtractRatio FUNCTION ---

   function [cmfRawTrace, cvfRegionTrace, cfRegionResponse, cnFramesInSample, cvfPixelResponse] = ...
         fhExtractRatioFun(fsData, cvnPixels, vnFrames, vnChannels, bUsedRR)
      
      if (~iscell(cvnPixels))
         cvnPixels = {cvnPixels};
      end
      
      nNumROIs = numel(cvnPixels);

      % - Concatenate pixels to extract
      cvnPixels = cellfun(@(c)(reshape(c, 1, [])), cvnPixels, 'UniformOutput', false);
      vnROISizes = cellfun(@numel, cvnPixels);
      mnROIBoundaries = [1 cumsum(vnROISizes)];
      mnROIBoundaries = [mnROIBoundaries(1:end-1)' mnROIBoundaries(2:end)'];
      vnExtractPixels = [cvnPixels{:}];

      % - Extract data from stack
      tfTrace = double(fsData(vnExtractPixels, vnFrames, vnChannels));
      
      % - Filter trace to remove zeros
      tfTrace(tfTrace == 0) = nan;
      
      % - Compute ratio
      mfRawTrace = tfTrace(:, :, 1) ./ tfTrace(:, :, 2);
      
      % - Calculate deltaR/R
      if (bUsedRR)
         mfBlankTrace = double(fsData.BlankFrames(vnExtractPixels, vnFrames));
         mfRawTraceDRR = (mfRawTrace - mfBlankTrace) ./ mfBlankTrace;
         mfRawTraceDRR(isnan(mfBlankTrace)) = mfRawTrace(isnan(mfBlankTrace));
         mfRawTrace = mfRawTraceDRR;
      end
      
      % - Extract regions
      for (nROI = nNumROIs:-1:1)
         vnThesePixels = mnROIBoundaries(nROI, 1):mnROIBoundaries(nROI, 2);
         
         cmfRawTrace{nROI} = mfRawTrace(vnThesePixels, :);
         
         if (nargout > 1)
            cvfRegionTrace{nROI} = nanmean(cmfRawTrace{nROI}, 1);
         end
         
         if (nargout > 2)
            cfRegionResponse{nROI} = nanmean(cvfRegionTrace{nROI});
         end
         
         if (nargout > 3)
            cnFramesInSample{nROI} = numel(vnFrames);
         end
         
         if (nargout > 4)
            cvfPixelResponse{nROI} = nanmean(cmfRawTrace{nROI}, 2);
         end
      end
      
      if (nNumROIs == 1)
         cmfRawTrace = cmfRawTrace{1};
         
         if (nargout > 1)
            cvfRegionTrace = cvfRegionTrace{1};
         end
         
         if (nargout > 2)
            cfRegionResponse = cfRegionResponse{1};
         end
         
         if (nargout > 3)
            cnFramesInSample = cnFramesInSample{1};
         end
         
         if (nargout > 4)
            cvfPixelResponse = cvfPixelResponse{1};
         end
      end      
   end
end
