function [mfFrameOffsets] = Align(oStack, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM)

% Align - METHOD Calculate an alignment for a stack
%
% Usage: [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM>)
%        [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling, nReferenceFrame, nWindowLength, vfSpatFreqCutoffCPUM>)
%
% The stack 'oStack' will be aligned, and the alignment will be recorded in the
% stack.

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

% - Compute frame alignments
mfFrameOffsets = GetStackAlignment( oStack, vnChannel, bProgressive, nUpsampling, ...
                                    mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM);
oStack.mfFrameShifts = mfFrameOffsets;

% - Clear aligned frame cache
oStack.vbCachedAlignedFrames = [];
oStack.oAlignedFrameCache = [];

% - Restore dF/F convertion
oStack.bConvertToDFF = bConvertToDFF;

% --- END of Align METHOD ---
