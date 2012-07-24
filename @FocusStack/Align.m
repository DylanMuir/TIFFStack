function [mfFrameOffsets] = Align(oStack, varargin)

% Align - METHOD Calculate an alignment for a stack
%
% Usage: [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength>)
%        [<mfFrameOffsets>] = Align(oStack <, vnChannel, bProgressive, nUpsampling, nReferenceFrame, nWindowLength>)
%
% The stack 'oStack' will be aligned, and the alignment will be recorded in the
% stack.

% -- Turn off dF/F conversion
bConvertToDFF = oStack.bConvertToDFF;
oStack.bConvertToDFF = false;

% - Was a reference supplied?
if (nargin > 4)
   % - Was it a reference frame index?
   if (isscalar(varargin{4}))
      % - Yes, so extract the requisite frame (mean of all channels)
      varargin{4} = nanmean(oStack.ExtractFrames({':', ':', varargin{4}, varargin{1}}), 4);
   end
end

% - Compute frame alignments
mfFrameOffsets = GetStackAlignment(oStack, varargin{:});
oStack.mfFrameShifts = mfFrameOffsets;

% - Clear aligned frame cache
oStack.vbCachedAlignedFrames = [];
oStack.oAlignedFrameCache = [];

% - Restore dF/F convertion
oStack.bConvertToDFF = bConvertToDFF;

% --- END of Align METHOD ---
