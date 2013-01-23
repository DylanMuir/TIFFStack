function [mfFrameOffsets] = Align(oStack, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM, fMaxShiftPerFrame)

% Align - METHOD Calculate an alignment for a stack
%
% Usage: [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM, fMaxShiftPerFrame>)
%        [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling, nReferenceFrame, nWindowLength, vfSpatFreqCutoffCPUM, fMaxShiftPerFrame>)
%        [<mfFrameOffsets>] = Align(oStack <, cvnChannelComb, ...>)
%
% The stack 'oStack' will be aligned, and the alignment will be recorded in the
% stack.
%
% The optional argument 'vnChannel' specifies which stack channel to use for
% registration.  If you want to use several channels summed, specify a vector of
% channel indices.
%
% 'cvnChannelComb' can be provided instead; this is a cell array, where the
% first cell contains a function handle that will be applied to a tensor [X Y F
% C], and must produce an output [X Y F] by combining channels in some way.  The
% second cell contains a vector of channel indices to use.  For example:
%    {@(t)nansum(t, 4) [1 2]}
% would extract channels 1 and 2 then sum them, and perform the alignment on
% the resulting image.
%
% The optional argument 'bProgressive' determines whether shifts should be
% calculated between successive frames ('bProgressive' = true), or whether all
% shifts should be calculated with the first frame ('bProgressive' = false,
% default).
%
% The optional argument 'nUpsampling' specifies that the registration should be
% computed to 1/'nUpsampling' pixels.  Default is 1, meaning that registration
% occurs to single pixel resolution.
%
% 'mfReferenceImage' is an optional reference image used for alignment, rather
% than the initial stack frame.  Optionally, a scalar frame index can be
% supplied instead.  In this case, the indicated frame will be extracted from
% the stack to be used as a reference for alignment.
%
% 'nWindowLength' is an optional parameter that specifies the number of frames
% to average together, in a moving window, to determine the alignment for the
% current frame.  Default is 1, meaning that only a single frame is used.
%
% 'vfSpatFreqCutoffCPUM' is an optional parameter that defines the cutoff spatial
% frequencies to include in computing the mis-alignment.  The vector is
% [fMinFreqCPUM fMaxFreqCPUM], both in cycles per um.  Spatial frequencies
% between these limits are included by a band-pass filter.
%

% -- Turn off dF/F conversion

bConvertToDFF = oStack.bConvertToDFF;
oStack.bConvertToDFF = false;

% -- Check arguments

if (nargin < 1)
   help FocusStack/Align;
   error('FocusStack:Usage', '*** FocusStack/Align: Incorrect usage.');
end

if (~exist('vnChannel', 'var'))
   vnChannel = [];
end

if (~exist('bProgressive', 'var'))
   bProgressive = [];
end

if (~exist('nUpsampling', 'var'))
   nUpsampling = [];
end

if (~exist('mfReferenceImage', 'var'))
   mfReferenceImage = [];
end

if (~exist('nWindowLength', 'var'))
   nWindowLength = [];
end

if (~exist('vfSpatFreqCutoffCPUM', 'var'))
   vfSpatFreqCutoffCPUM = [];

elseif (~isempty(oStack.fPixelsPerUM))
   % - Modify cutoff frequencies to per-pixel frequencies
   vfSpatFreqCutoffCPUM = vfSpatFreqCutoffCPUM ./ oStack.fPixelsPerUM;
end

if (~exist('fMaxShiftPerFrame', 'var'))
   fMaxShiftPerFrame = inf;
end

% - Find trial starts and ends
w = warning('off', 'FocusStack:IncompleteInformation');
[nul, vnTrialIndex] = oStack.FrameStimulusInfo;
nNumTrials = numel(oStack.cstrFilenames);
mnTrialRanges = nan(nNumTrials, 2);

for (nTrialIndex = 1:nNumTrials)
   mnTrialRanges(nTrialIndex, :) = [find(vnTrialIndex == nTrialIndex, 1, 'first') find(vnTrialIndex == nTrialIndex, 1, 'last')];
end

% - Clear aligned frame cache
oStack.vbCachedAlignedFrames = [];
oStack.oAlignedFrameCache = [];
oStack.mfFrameShifts = zeros(size(oStack, 3), 2);

% - Compute frame alignments
mfFrameOffsets = GetStackAlignment( oStack, vnChannel, bProgressive, nUpsampling, ...
                                    mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM, ...
                                    mnTrialRanges);

%% -- Filter frame shifts

vfLastDists = [];
vfDists = nan;

while (~isequal(vfLastDists, vfDists))
   vfLastDists = vfDists;
   
   % - Compute Euclidean shifts
   vfDists = sqrt(sum(diff(mfFrameOffsets, 1, 1).^2, 2));
   
   % - Ignore the beginning of each trial
   vfDists(mnTrialRanges(:)) = 0;

   % - Which frame shifts should we ignore?
   vbIgnoreShift = vfDists > fMaxShiftPerFrame;
   
   while (any(vbIgnoreShift))
      % - Rplace the first found frame
      nReplaceInd = find(vbIgnoreShift, 1, 'first');
      mfFrameOffsets(nReplaceInd+1, :) = mfFrameOffsets(nReplaceInd, :);
      
      % - Fix up frame shift distances
      vfDists(nReplaceInd) = 0;
      if (nReplaceInd < numel(vfDists))
         vfDists(nReplaceInd+1) = sqrt(sum((mfFrameOffsets(nReplaceInd+2, :) - mfFrameOffsets(nReplaceInd+1, :)).^2, 2));
      end
      
      plot(vfDists);drawnow;

      % - Which frame shifts should we ignore?
      vbIgnoreShift = vfDists > fMaxShiftPerFrame;
   end
end
      
% - Assign frame shifts
oStack.mfFrameShifts = mfFrameOffsets;


% - Restore dF/F convertion
oStack.bConvertToDFF = bConvertToDFF;
warning(w);

% --- END of Align METHOD ---
