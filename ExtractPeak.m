function [fhExtractionFun] = ExtractPeak(nSampleFrames, nChannel, bUsedFF)

% ExtractPeak - Extraction function for returning the peak of the response
%
% Usage: [fhExtractionFun] = ExtractPeak(<nSampleFrames, nChannel, bUsedFF>)
% Usage: [fhExtractionFun] = ExtractPeak(<vnSampleOffsets, nChannel, bUsedFF>)
%
% This function will extract the trace of a region, find the peak point
% during that trace, and return the average of a few frames around that
% point.
%
% 'nSampleFrames' specifies the number of frames to sample for the mean of
% the peak (default: 5).  This will correspond to a sampling window of
% [-2 -1 0 1 2] frames around the peak (placed at 0).
%
% Optionally, 'vnSampleOffsets' can be provided as a vector of offsets from
% the peak response location.  For example, [-2 -1 0 1 2] would result in
% the default sampling window.  [0 1 2] would sample the peak and the two
% following frames.
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

if (numel(nSampleFrames) == 1)
   vnSampleWindow = (-floor(nSampleFrames/2)):(ceil(nSampleFrames/2)-1);
else
   vnSampleWindow = nSampleFrames;
end


% -- Return function handle

fhExtractionFun = @(fsData, vnPixels, vnFrames)fhExtractPeak(fsData, vnPixels, vnFrames, vnSampleWindow, nChannel, bUsedFF);


% --- END of ExtractPeak.m ---

   function [cmfRawTrace, cvfRegionTrace, cfRegionResponse, cnFramesInSample, cvfPixelResponse] = ...
         fhExtractPeak(fsData, cvnPixels, vnFrames, vnSampleWindow, nChannel, bUsedFF)
      
      if (~iscell(cvnPixels))
         cvnPixels = {cvnPixels};
      end
      
      nNumROIs = numel(cvnPixels);

      % - Convert logical indexing to numerical indexing
      vbIsLogical = cellfun(@islogical, cvnPixels);
      cvnPixels(vbIsLogical) = cellfun(@(c)(find(c)), cvnPixels(vbIsLogical), 'UniformOutput', false);

      if (islogical(vnFrames))
         vnFrames = find(vnFrames);
      end
      
      % - Concatenate pixels to extract
      cvnPixels = cellfun(@(c)(reshape(c, 1, [])), cvnPixels, 'UniformOutput', false);
      vnROISizes = cellfun(@numel, cvnPixels);
      mnROIBoundaries = [1 cumsum(vnROISizes)];
      mnROIBoundaries = [mnROIBoundaries(1:end-1)' mnROIBoundaries(2:end)'];
      vnExtractPixels = [cvnPixels{:}];

       % - Extract data from stack
      mfRawTrace = double(fsData(vnExtractPixels, vnFrames, nChannel));
      if (bUsedFF)
         mfBlankTrace = double(fsData.BlankFrames(vnExtractPixels, vnFrames));
      end
      
      % - Extract regions
      for (nROI = nNumROIs:-1:1)
         vnThesePixels = mnROIBoundaries(nROI, 1):mnROIBoundaries(nROI, 2);
         vfThisRawTrace = nanmean(mfRawTrace(vnThesePixels, :), 1);
         
         % - Calculate deltaF/F
         if (bUsedFF)
            vfThisBlankTrace = nanmean(mfBlankTrace(vnThesePixels, :), 1);
            vfThisRawTraceDFF = (vfThisRawTrace - vfThisBlankTrace) ./ vfThisBlankTrace;
            vfThisRawTraceDFF(isnan(vfThisBlankTrace)) = vfThisRawTrace(isnan(vfThisBlankTrace));
            vfThisRawTrace = vfThisRawTraceDFF;
         end
         
         cmfRawTrace{nROI} = mfRawTrace(vnThesePixels, :);
         
         if (nargout > 1)
            cvfRegionTrace{nROI} = vfThisRawTrace;
         end
         
         if (nargout > 2)
            % - Locate peak
            [nul, nPeakIndex] = nanmax(cvfRegionTrace{nROI});
            
            % - Average frames surrounding peak, according to sampling window
            vnWindow = nPeakIndex + vnSampleWindow;
            vnWindow = vnWindow(vnWindow >= 1);
            vnWindow = vnWindow(vnWindow <= numel(cvfRegionTrace{nROI}));
            
            cfRegionResponse{nROI} = nanmean(cvfRegionTrace{nROI}(vnWindow));
         end
         
         if (nargout > 3)
            cnFramesInSample{nROI} = numel(vnWindow);
         end
         
         if (nargout > 4)
            cvfPixelResponse{nROI} = nanmax(cmfRawTrace{nROI}, [], 2);
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