function [vtGlobalTime, ...
          vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
          vnStimulusSeqID, vtTimeInStimPresentation, ...
          vnPresentationIndex, vbUseFrame] = ...
            FrameStimulusInfo(oStack, vnFrameIndices)

% FrameStimulusInfo - METHOD Return stimulus information about stack frames
%
% Usage: [vtGlobalTime, ...
%         vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
%         vnStimulusSeqID, vtTimeInStimPresentation, ...
%         vnPresentationIndex, vbUseFrame] = ...
%           FrameStimulusInfo(oStack <, vnFrameIndices>)
%
% 'oStack' is a FocusStack.
%
% 'vnFrameIndices' is a vector of frame indices into 'oStack', for which we
% would like some extra information.  If not provided, information for all
% frames in order will be returned.
%
% 'vtGlobalTime' will be a vector of time, in seconds, corresponding to each
% frame in 'vnFrameIndices'.
%
% 'vnBlockIndex' will be a vector of block IDs, corresponding to each frame in
% 'vnFrameIndices'.
%
% 'vnFrameInBlock' will be a vector of frame indices into the enclosing block,
% for each frame in 'vnFrameIndices'.
%
% 'vtTimeInBlock' will be a vector of time points, in seconds, as an offset from
% the start of each enclosing block, for each frame in 'vnFrameIndices'.
%
% 'vnStimulusSeqID' will be a vector of stimulus sequence IDs corresponding to
% each frame in 'vnFrameIndices'.
%
% 'vtTimeInStimPresentation' will be a vector of time offsets, in seconds, from
% the start of the enclosing stimulus presentation.
%
% 'vnPresentationIndex' will be a vector of indices, from 1 to N (where N is the
% total number of stimulus presentations), corresponding to each frame in
% 'vnFrameIndices'.
%
% 'vbUseFrame' will be a vector of booleans, indicating whether the
% corresponding frame in 'vnFrameIndices' should be used for analysis, according
% to the 'mtStimulusUseTimes' data for the stack.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 13th January, 2011

% -- Check for valid indices

if (nargin < 1)
   error('FocusStack:InvalidUsage', ...
      '*** FocusStack/FrameStimulusInfo: Invalid usage.');
end

nNumFrames = size(oStack, 3);

if (~exist('vnFrameIndices', 'var') || isempty(vnFrameIndices))
   vnFrameIndices = 1:nNumFrames;

elseif (any(vnFrameIndices(:) < 1) || any(vnFrameIndices(:) > nNumFrames))
   error('FocusStack:InvalidArgument', ...
      '*** FocusStack/FrameStimulusInfo: ''vnFrameIndices'' must be limited to [1 .. %d] for this stack.', ...
      nNumFrames);
end


% -- Does the block have the required information?

if (isempty(oStack.tFrameDuration))
   error('FocusStack:IncompleteInformation', ...
      '*** FocusStack/FrameStimulusInfo: The frame rate was not available for this stack.');
end

if (isempty(oStack.vtStimulusStartTimes) || isempty(oStack.vtStimulusDurations))
   warning('FocusStack:IncompleteInformation', ...
      '--- FocusStack/FrameStimulusInfo: Not all the required information is present in the stack.');
   bComputeStimInfo = false;
else
   bComputeStimInfo = true;
end


% -- Calculate global time

vtGlobalTime = (vnFrameIndices-1) * oStack.tFrameDuration;


% -- Work out which block(s) we're in, which frame in which block

nNumBlocks = numel(oStack.cstrFilenames);
% nBlockLength = nNumFrames / nNumBlocks;
vnBlockFrameIndices = vnFrameIndices;

for (nBlock = 1:nNumBlocks)
   vbInBlock = vnFrameIndices <= sum(oStack.vnNumFrames(1:nBlock));
   vbInBlock = vbInBlock & (vnFrameIndices > sum(oStack.vnNumFrames(1:(nBlock-1))));
   vnBlockIndex(vbInBlock) = nBlock; %#ok<AGROW>
   vnFrameInBlock(vbInBlock) = mod(vnFrameIndices(vbInBlock)-1, oStack.vnNumFrames(nBlock))+1; %#ok<AGROW>
   vnBlockFrameIndices(vbInBlock) = nan;
end

vtTimeInBlock = (vnFrameInBlock-1) * oStack.tFrameDuration;


% -- Compute stimulus information

if (bComputeStimInfo)
   % - Get stimulus start times and nominal durations
   vtStimulusStartTimes = oStack.vtStimulusStartTimes;
   vtStimulusEndTimes = oStack.vtStimulusEndTimes;
   vtStimulusDurations = oStack.vtStimulusDurations;
   mtStimulusUseTimes = oStack.mtStimulusUseTimes;
   
   % - Compute time in stimulus segment
   cvnSequenceIDs = oStack.cvnSequenceIDs;
   vnStimOrder = vertcat(cvnSequenceIDs{:});
   nNumPresentations = numel(vnStimOrder);
   vnStimulusSeqID = nan(1, numel(vnFrameIndices));
   vnPresentationIndex = nan(1, numel(vnFrameIndices));
   vtTimeInStimPresentation = nan(1, numel(vnFrameIndices));
   vbUseFrame = false(1, numel(vnFrameIndices));

   vbAssignedStimID = false(1, numel(vnFrameIndices));
   nStimPresNum = 1;
   while (any(~vbAssignedStimID) && (nStimPresNum <= nNumPresentations))
      vbInsideThisStimulus =  (vtGlobalTime >= vtStimulusStartTimes(nStimPresNum)) & ...
                              (vtGlobalTime <= vtStimulusEndTimes(nStimPresNum)); 
      vbAssignedStimID = vbAssignedStimID | vbInsideThisStimulus;
      
      vnStimulusSeqID(vbInsideThisStimulus) = vnStimOrder(nStimPresNum);
      vtTimeInStimPresentation(vbInsideThisStimulus) = abs(vtGlobalTime(vbInsideThisStimulus) - vtStimulusStartTimes(nStimPresNum));
      vnPresentationIndex(vbInsideThisStimulus) = nStimPresNum;
      
      % - Which frames should be used for analysis?
      if (~isnan(vnStimOrder(nStimPresNum)))
         vtThisUseTimes = mtStimulusUseTimes(vnStimOrder(nStimPresNum), :);
         vbUseFrame(vbInsideThisStimulus) =  (vtTimeInStimPresentation(vbInsideThisStimulus) >= vtThisUseTimes(1)) & ...
                                             (vtTimeInStimPresentation(vbInsideThisStimulus) <= vtThisUseTimes(2));
      end
      
      % - Move to next stimulus presentation
      nStimPresNum = nStimPresNum + 1;
   end
   
else
   % - Return empty matrices for stimulus information
   vnStimulusSeqID = [];
   vtTimeInStimPresentation = [];
   vnPresentationIndex = [];
   vbUseFrame = [];
   
end

% --- END of FrameStimulusInfo METHOD ---
