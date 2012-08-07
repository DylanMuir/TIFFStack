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

   function [cmfRawTrace, cvfRegionTrace, cfRegionResponse, cnFramesInSample, cvfPixelResponse] = ...
         fhExtractMeanFun(fsData, cvnPixels, vnFrames, nChannel, bUsedFF)
      
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
      
%       if (numel(find(vnFrames)) == size(fsData, 3))
%          for (nFrame = numel(vnFrames):-1:1)
%             mfRawTrace(:, nFrame, 1) = double(fsData(vnPixels, vnFrames(nFrame), nChannel));
%          end
%       else
         mfRawTrace = double(fsData(vnExtractPixels, vnFrames, nChannel));
%       end
      
      % - Calculate deltaF/F
      if (bUsedFF)
         mfBlankTrace = double(fsData.BlankFrames(vnExtractPixels, vnFrames));
         mfRawTraceDFF = (mfRawTrace - mfBlankTrace) ./ mfBlankTrace;
         mfRawTraceDFF(isnan(mfBlankTrace)) = mfRawTrace(isnan(mfBlankTrace));
         mfRawTrace = mfRawTraceDFF;
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

% --- END of ExtractMean.m ---
