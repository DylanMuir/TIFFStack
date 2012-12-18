function [vhFigures, mfRegionResponses] = ...
            PlotRegionResponses(  mfRegionResponses, fsStack, sRegions, ...
                                  mtStimLabelTimes, ...
                                  fDFFScaleBarLength, tTimeScaleBarLength, ...
                                  nLimitNumberOfCells, cstrAnnotation, ...
                                  mbHighlightStims)

% PlotRegionResponses - FUNCTION Make a set of plots showing the response traces of identified regions
%
% Usage: [vhFigures, mfRegionResponses] = ...
%           PlotRegionResponses(  <mfRegionResponses>, fsStack, sRegions, ...
%                                 <mtStimLabelTimes>, ...
%                                 <fDFFScaleBarLength, tTimeScaleBarLength, ...
%                                 nLimitNumberOfCells, cstrAnnotation, ...
%                                 mbHighlightStims>)
%
% If 'fsStack' is supplied, region traces will be extracted and returned in
% 'mfRegionResponses'.  The traces for channel 1 only will be used.  The region
% traces can be supplied as an argument in 'mfRegionResponses', in which case
% they should be in a 2D NxT matrix, where N is the number of regions in
% 'sRegions' and T is the full time length of the source stack.
%
% 'sRegions' is a regions structure as extracted by the bwconncomp function.
%
% 'mtStimLabelTimes' is an optional parameter that specifies times to label for
% each stimulus, if desired.  This is an Sx2 array.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 7th December, 2010

% -- Defaults

DEF_nNumResponsesPerFigure = 4;


% -- Check arguments

if (nargin < 3)
   disp('*** PlotRegionResponses: Incorrect usage');
   help PlotRegionResponses;
   vhFigures = [];
   return;
end

if (~exist('mbHighlightStims', 'var') || isempty(mbHighlightStims))
   mbHighlightStims = false(size(mfRegionResponses));
end

% - Check that required stack data are present
if (  isempty(fsStack.cvnSequenceIDs) || ...
      isempty(fsStack.vtStimulusDurations) || ...
      isempty(fsStack.tFrameDuration))
   disp('--- PlotRegionResponses: Warning: Not all required stack data is available.');
   bSegmentStack = false;
else
   bSegmentStack = true;
end

% - Check if we need to extract region responses
if (isempty(mfRegionResponses))
   % - We should extract region responses and average
   disp('--- PlotRegionResponses: Extracting region traces...');

   % - Extract activity traces for all regions (assume channel 1)
   mfPixelTraces = fsStack.AlignedStack(vertcat(sRegions.PixelIdxList{:}), :, 1);
   
   mfRegionResponses = zeros(sRegions.NumObjects, size(fsStack, 3));
   vnNumTracesRegion = cellfun(@numel, sRegions.PixelIdxList);
   
   for (nRegion = 1:sRegions.NumObjects)
      % - Extract the region traces
      nFirstTrace = sum(vnNumTracesRegion(1:nRegion-1))+1;
      nLastTrace = sum(vnNumTracesRegion(1:nRegion));
      
      % - Average together all traces for this region
      mfRegionResponses(nRegion, :) = mean(mfPixelTraces(nFirstTrace:nLastTrace, :), 1);
   end
end

% - Get stack information and check
nNumFrames = size(fsStack, 3);
cvnSequenceIDs = fsStack.cvnSequenceIDs;
vnStimOrder = vertcat(cvnSequenceIDs{:});
nNumStimuli = fsStack.nNumStimuli;
nNumRegions = sRegions.NumObjects;
nNumBlocks = numel(fsStack.cstrFilenames);

if (~isequal(size(mfRegionResponses, 1), nNumRegions))
   error('PlotRegionResponses:BadArguments', ...
         '*** PlotRegionResponses: ''mfRegionResponses'' must contain a row for each region in ''sRegions''.');
end

if (~isequal(size(mfRegionResponses, 2), nNumFrames))
   error('PlotRegionResponses:BadArguments', ...
         '*** PlotRegionResponses: ''mfRegionResponses'' must contain a column for each frame in ''fsStack''.');
end

if (~exist('mtStimLabelTimes', 'var'))
   mtStimLabelTimes = [];
end

if (~isempty(mtStimLabelTimes) && ~isequal(size(mtStimLabelTimes, 1), nNumStimuli))
   error('PlotRegionResponses:BadArguments', ...
         '*** PlotRegionResponses: ''mtStimLabelTimes'' must have the correct number of stimuli.');
end


% -- Extract stimulus information for each frame

[vtGlobalTime, ...
 vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
 vnStimulusSeqID, vtTimeInStimPresentation] = ...
 	FrameStimulusInfo(fsStack, 1:nNumFrames);


% -- Decide where to place stimuli in "time"

if (bSegmentStack)
   vtStimDurations = fsStack.vtStimulusDurations;
   vtStimulusBlockStartTimes = cumsum([0; vtStimDurations]);
else
   vtStimulusBlockStartTimes = 0;
   nNumStimuli = 1;
   vnStimulusSeqID = ones(size(vtGlobalTime));
   vtTimeInStimPresentation = repmat(vtGlobalTime(vnBlockIndex == 1), 1, nNumBlocks);
end


% -- Make some plots

vhFigures = [];

if (exist('nLimitNumberOfCells', 'var') && ~isempty(nLimitNumberOfCells))
   if (isscalar(nLimitNumberOfCells))
      vnRegions = 1:min(nNumRegions, nLimitNumberOfCells);
      
   else
      vnRegions = nLimitNumberOfCells;
   end
end

% - Loop over regions
nNumRegions = numel(vnRegions);
for (nRegion = 1:nNumRegions)
   nRegionToPlot = vnRegions(nRegion);
   if (mod(nRegion-1, DEF_nNumResponsesPerFigure) == 0)
      vhFigures(end+1) = figure; %#ok<AGROW>
      set(gcf, 'Color', 'w');
   end
   
   % - Create a subplot
   subplot(DEF_nNumResponsesPerFigure/2, DEF_nNumResponsesPerFigure/2, mod(nRegion-1, DEF_nNumResponsesPerFigure)+1);
   
   % - Determine plot limits
   fYMin = min(mfRegionResponses(nRegionToPlot, :));
   fYMax = max(mfRegionResponses(nRegionToPlot, :));
   vfYLims = [fYMin fYMax];

   % - Loop over stimulus segments
   for (nStim = 1:nNumStimuli)
      % - Label the stimulus presentation time, if requested
      if (~isempty(mtStimLabelTimes) && ~any(isnan(mtStimLabelTimes(nStim, :))))
         tStimStart = vtStimulusBlockStartTimes(nStim) + mtStimLabelTimes(nStim, 1);
         tStimEnd = vtStimulusBlockStartTimes(nStim) + mtStimLabelTimes(nStim, 2);
         if size(mtStimLabelTimes, 2) == 3;
             tAdaptationStart = vtStimulusBlockStartTimes(nStim) + mtStimLabelTimes(nStim, 3);
         else
             tAdaptationStart = [];
         end
         vtStimTimes = [tStimStart tStimEnd tAdaptationStart];
         
         if (mbHighlightStims(nRegion, nStim))
            vfColour = 0.75 * [1 0 0];
         else
            vfColour = 0.8 * [1 1 1];
         end
         
         hPatchStim = fill(vtStimTimes([1 1 2 2 1]), vfYLims([1 2 2 1 1]), vfColour);
         set(hPatchStim, 'LineStyle', 'none');
         hold on;
         if size(mtStimLabelTimes, 2) == 3;
             hPatchAdapt = fill(vtStimTimes([3 3 1 1 3]), vfYLims([1 2 2 1 1]), ([1 1 1] + vfColour) / 2);
             set(hPatchAdapt, 'LineStyle', 'none');
         end
         
      end      
      
      % - Plot individual trials for this stimulus
      vfAvgStimTrace = [];
      vnNorm = [];
      for (nBlock = 1:nNumBlocks)
         % - Find frames corresponding to this stimulus ID in this block 
         vbStimBlockFrames = (vnStimulusSeqID == nStim) & (vnBlockIndex == nBlock);
         vfThisTrialTrace = mfRegionResponses(nRegionToPlot, vbStimBlockFrames);
         
         % - Find times to plot this response trace on
         vtThisPresTimes = vtStimulusBlockStartTimes(nStim) + vtTimeInStimPresentation(vbStimBlockFrames);

         if (any(isnan(vfThisTrialTrace)))
            disp('bluh!');
         end
         
         % - Plot this trial trace
         plot(vtThisPresTimes, vfThisTrialTrace, 'Color', 0.6 * [1 1 1], 'LineWidth', 1);
         hold on;
            
         % - Accumulate stim traces
         nStimLength = numel(vfThisTrialTrace);
         if (nStimLength > numel(vfAvgStimTrace))
            vfAvgStimTrace(nStimLength) = 0; %#ok<AGROW>
            vnNorm(nStimLength) = 0; %#ok<AGROW>
         end
         
         vfAvgStimTrace(1:nStimLength) = vfAvgStimTrace(1:nStimLength) + vfThisTrialTrace(1:nStimLength);
         vnNorm(1:nStimLength) = vnNorm(1:nStimLength) + 1;
      end
      
      if (any(isnan(vfAvgStimTrace ./ vnNorm)))
         disp('bluh!');
      end
      
      % - Plot the average fluorescence trace for the region
      vbPlotTimePoints = vnNorm == max(vnNorm);
      nStimLength = numel(vfAvgStimTrace);
      vtThisStimTime = vtStimulusBlockStartTimes(nStim) + (0:nStimLength-1) * fsStack.tFrameDuration;
      plot(vtThisStimTime(vbPlotTimePoints), vfAvgStimTrace(vbPlotTimePoints)...
          ./ vnNorm(vbPlotTimePoints), 'Color', 0 * [1 1 1], 'LineWidth', 2);
   end
   
   % - Get axis limits
   axis tight;
   vfLims = axis;
   
   if (any(isnan(vfYLims)) || (vfYLims(1) == vfYLims(2)))
      vfYLims = vfLims(3:4);
   end
   
   % - Plot a dF/F scale bar
   tBarPos = vfLims(1)+5;
   if (exist('fDFFScaleBarLength', 'var'))
      vfBarLims = 0.8*max(vfYLims) + [-1 0] * fDFFScaleBarLength;
   else
      vfBarLims = 0.8*max(vfYLims) + [-1 0] * 0.1;
   end
   plot([tBarPos tBarPos], vfBarLims, 'k-', 'LineWidth', 5);
   
   % - Plot a time scale bar
   if (exist('tTimeScaleBarLength', 'var'))
      vfTBarLims = tBarPos + [0 tTimeScaleBarLength];
      plot(vfTBarLims, vfBarLims([2 2]), 'k-', 'LineWidth', 5);
   end
   
   % - Set the plot properties
   axis([0 vfLims(2) vfYLims]);
%    axis tight;
   box off;
   
   % plot only every "DEF_nNumResponsesPerFigure"-th figure
%    if (nRegion ~= sRegions.NumObjects) && (mod(nRegion-1, DEF_nNumResponsesPerFigure) ~= (DEF_nNumResponsesPerFigure-1))
%     set(gca, 'XTick', []);
%    end
   
   set(gca, 'Color', 'none', 'YTick', [], 'LineWidth', 1, 'FontSize', 18);
   ylabel(gca, nRegionToPlot, 'FontSize', 24);
   
   % - Add the annotation, if required
   if (exist('cstrAnnotation', 'var') && ~isempty(cstrAnnotation))
      text(vfLims(2), 0, cstrAnnotation{nRegion}, 'FontSize', 12);
   end
end



% --- END of PlotRegionResponses.m ---
