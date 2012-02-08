function [cvnSegments, cvnResponseSegments, cvnBlockEndSegments] = ...
   SegmentStack(tFrameRate, vtStimSegmentDurations, vnStimOrder, mtUseStimDurations, nFramesPerBlock)

% SegementStack - FUNCTION Compute stack frame indicies for a set of stimuli
%
% Usage: [cvnSegemnts, cvnResponseSegments, cvnBlockEndSegments] = ...
%           SegmentStack(tFrameRate, vtStimSegmentDurations, vnStimOrder <, mtUseStimDurations, nFramesPerBlock>)
%
% A value of NaN in 'vnStimOrder' implies that the rest of the block (where each
% block has 'nFramesPerBlock' frames) will be skipped.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 26th November, 2010

% -- Check arguments

if (nargin < 3)
   disp('*** SegmentStack: Incorrect usage');
   help SegmentStack;
   return;
end

nNumStimuli = numel(vtStimSegmentDurations);

% - Use all stimulus frames, if not otherwise specified
if (~exist('mtUseStimDurations', 'var') || isempty(mtUseStimDurations))
   mtUseStimDurations(:, 1) = zeros(nNumStimuli, 1);
   mtUseStimDurations(:, 2) = reshape(vtStimSegmentDurations, [], 1);
end

% - Check the size of 'mtUseStimDurations'
if (size(mtUseStimDurations, 1) ~= nNumStimuli)
   disp('*** SegementStack: ''mtUseStimDurations'' must have the same number of rows as ''vtStimSegmentDurations''.');
   return;
end

% - Check whether "skip" stimuli were specified, and if so whether the block
% size was also specified
if (any(isnan(vnStimOrder(:))))
   if (~exist('nFramesPerBlock', 'var') || isempty(nFramesPerBlock))
      disp('*** ''nFramesPerBlock'' must be specifed if "skip" stimuli are specified.');
      return;
   end
end


% -- Carve up stimuli

% - Pre-allocate stimulus segments
cvnSegments = cell(nNumStimuli, 1);
cvnResponseSegments = cvnSegments;
cvnBlockEndSegments = {};

nNumPresentations = numel(vnStimOrder);
tTime = 0;


% for (nStimPresentation = 1:nNumPresentations)
vtSegmentStartTimes = vtStimSegmentDurations(vnStimOrder);



% - Compute time ranges for each stimulus presentation
for (nStim = 1:nNumPresentations)
   % - Check for a "skip to end of block" stimulus
   if (isnan(vnStimOrder(nStim)))
      tTimeLeftInBlock = (nFramesPerBlock/tFrameRate) - mod(tTime, nFramesPerBlock/tFrameRate);
      cvnBlockEndSegments{end+1} = round((tTime*tFrameRate+1)):((tTime + tTimeLeftInBlock)*tFrameRate); %#ok<AGROW>
      % - Skip to the next multiple of nFramesPerBlock/tFrameRate
      tTime = tTime + tTimeLeftInBlock;
      
   else
      % - This is a normal stimulus
      tThisStimDuration = vtStimSegmentDurations(vnStimOrder(nStim));
      
      % - Compute time and frame ranges for this stimulus
      vnThisStimFrames = floor((tTime*tFrameRate+1)):ceil(((tTime+tThisStimDuration)*tFrameRate));
      vnUseFrames = floor((mtUseStimDurations(vnStimOrder(nStim), 1)*tFrameRate+1)):ceil((mtUseStimDurations(vnStimOrder(nStim), 2)*tFrameRate));
      
      % - Record the frames for this stimulus
      cvnSegments{vnStimOrder(nStim)}{end+1}.vnFrames = vnThisStimFrames;
      cvnSegments{vnStimOrder(nStim)}{end}.vtTimeOffsets = vnThisStimFrames * tFrameRate
      
      cvnResponseSegments{vnStimOrder(nStim)}{end+1} = vnThisStimFrames(vnUseFrames);
      
      % - Advance time
      tTime = tTime + tThisStimDuration;
   end
end

% --- END of SegmentStack.m ---
