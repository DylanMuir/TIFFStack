function [mfFrameOffsets] = Align(oStack, varargin)

% Align - METHOD Calculate an alignment for a stack
%
% Usage: [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling>)
%
% The stack 'oStack' will be aligned, and the alignment will be recorded in the
% stack.

% -- Turn off dF/F conversion
bConvertToDFF = oStack.bConvertToDFF;
oStack.bConvertToDFF = false;

mfFrameOffsets = GetStackAlignment(oStack, varargin{:});
oStack.mfFrameShifts = mfFrameOffsets;
oStack.vbCachedAlignedFrames = [];
oStack.oAlignedFrameCache = [];

% - Restore dF/F convertion
oStack.bConvertToDFF = bConvertToDFF;

% --- END of Align METHOD ---
