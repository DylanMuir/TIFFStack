function [mfFrameOffsets] = AlignTrials(oStack, vnChannel, nUpsampling, mfReferenceImage, vfSpatFreqCutoffCPUM)

% AlignTrials - METHOD Correct mis-alignment on a trial-by-trial basis, not shifting each frame
%
% Usage: [<mfFrameOffsets>] = AlignTrials(oStack <, vnChannel, nUpsampling, mfReferenceImage, vfSpatFreqCutoffCPUM>)
%        [<mfFrameOffsets>] = AlignTrials(oStack <, vnChannel, nUpsampling, nReferenceFrame, vfSpatFreqCutoffCPUM>)
%
% Each trial in 'oStack' will be aligned as a block, and the estimated
% frame offsets returned in 'mfFrameOffsets'.  The alignment correction
% will also be applied to the stack.

% -- Defaults

DEF_vnChannel = 1;


% -- Turn off dF/F conversion

bConvertToDFF = oStack.bConvertToDFF;
oStack.bConvertToDFF = false;

% -- Check arguments

if (nargin < 1)
   help FocusStack/Align;
   error('FocusStack:Usage', '*** FocusStack/Align: Incorrect usage.');
end

if (~exist('vnChannel', 'var'))
   vnChannel = DEF_vnChannel;
end

if (~exist('nUpsampling', 'var'))
   nUpsampling = [];
end

% - Extract frame information
w = warning('off', 'FocusStack:IncompleteInformation');
[nul, vnBlockIndex] = oStack.FrameStimulusInfo;

if (~exist('mfReferenceImage', 'var'))
   % - Take the first trial as the reference
   mfReferenceImage = ExtractSummedFrames(oStack, {':', ':', vnBlockIndex == 1, vnChannel}, false) ./ nnz(vnBlockIndex == 1);
   nTrialStart = 2;
elseif (isscalar(mfReferenceImage))
   mfReferenceImage = oStack(:, :, mfReferenceImage, vnChannel);
   nTrialStart = 1;
end

if (~exist('vfSpatFreqCutoffCPUM', 'var'))
   vfSpatFreqCutoffCPUM = [];

elseif (~isempty(oStack.fPixelsPerUM))
   % - Modify cutoff frequencies to per-pixel frequencies
   vfSpatFreqCutoffCPUM = vfSpatFreqCutoffCPUM ./ oStack.fPixelsPerUM;
end

% - Clear aligned frame cache
oStack.vbCachedAlignedFrames = [];
oStack.oAlignedFrameCache = [];
oStack.mfFrameShifts = zeros(size(oStack, 3), 2);

% -- Loop over trials

mfFrameOffsets = zeros(size(oStack, 3), 2);

for (nTrialIndex = nTrialStart:numel(oStack.cstrFilenames))
   % - Average the frames for this trials
   vbThisBlock = vnBlockIndex == nTrialIndex;
   mfTrialAverage = ExtractSummedFrames(oStack, {':', ':', vbThisBlock, vnChannel}, false) ./ nnz(vbThisBlock);
   
   % - Estimate the misalignemnt for this trial
   mfFrameOffsets(vbThisBlock, :) = ...
      repmat(GetStackAlignment(mfTrialAverage, 1:size(mfTrialAverage, 3), false, nUpsampling, ...
                               mfReferenceImage, 1, vfSpatFreqCutoffCPUM), ...
             nnz(vbThisBlock), 1);
end

oStack.mfFrameShifts = mfFrameOffsets;


% - Restore dF/F conversion
oStack.bConvertToDFF = bConvertToDFF;
warning(w);


% --- END of AlignTrials.m ---
