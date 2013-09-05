function [mfFrameOffsets] = AlignTrials(oStack, vnChannel, nUpsampling, mfReferenceImage, vfSpatFreqCutoffCPUM)

% AlignTrials - METHOD Correct mis-alignment on a trial-by-trial basis, not shifting each frame
%
% Usage: [<mfFrameOffsets>] = AlignTrials(oStack <, vnChannel, nUpsampling, mfReferenceImage, vfSpatFreqCutoffCPUM>)
%        [<mfFrameOffsets>] = AlignTrials(oStack <, vnChannel, nUpsampling, nReferenceFrame, vfSpatFreqCutoffCPUM>)
%        [<mfFrameOffsets>] = AlignTrials(oStack <, vnChannel, nUpsampling, vnReferenceWindow, vfSpatFreqCutoffCPUM>)
%
% Each trial in 'oStack' will be aligned as a block, and the estimated
% frame offsets returned in 'mfFrameOffsets'.  The alignment correction
% will also be applied to the stack.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 2012

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
   % - Assume 'mfReferenceImage' is a frame reference
   mfReferenceImage = oStack(:, :, mfReferenceImage, vnChannel);
   nTrialStart = 1;
   
elseif (any(size(mfReferenceImage) == 1))
   % - Assume 'mfReferenceImage' is a list of frames to average for a reference image
   vnRefWindow = mfReferenceImage;
   mfReferenceImage = ExtractSummedFrames(oStack, {':', ':', vnRefWindow, vnChannel}, false) ./ numel(mfReferenceImage);
   nTrialStart = 1;
   
else
   % - Assume 'mfReferenceImage' is really an image to align to
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
   
   if (exist('vnRefWindow', 'var'))
      % - Use a referencing window
      vnThisBlock = find(vbThisBlock);
      vnThisWindow = vnThisBlock(vnRefWindow);
      mfTrialAverage = ExtractSummedFrames(oStack, {':', ':', vnThisWindow, vnChannel}, false) ./ nnz(vbThisBlock);
   else
      % - Take the whole trials
      mfTrialAverage = ExtractSummedFrames(oStack, {':', ':', vbThisBlock, vnChannel}, false) ./ nnz(vbThisBlock);
   end
   
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
