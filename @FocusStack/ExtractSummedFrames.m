function [tfData] = ExtractSummedFrames(oStack, sSubs, bAligned)

% ExtractSummedFrames - METHOD Extract summed frames from a stack

% -- Check whether the reference is even possible

if (iscell(sSubs))
   subs = sSubs;
   sSubs = [];
   sSubs.type = '()';
   sSubs.subs = subs;
end

if (numel(sSubs) > 1)
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractSummedFrames: Only one referencing level is possible.');
end

if (~isequal(sSubs.type, '()'))
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractSummedFrames: Only ''()'' referencing is supported.');
end


% -- Check arguments

if (~exist('bAligned', 'var') || isempty(bAligned))
   bAligned = false;
end


% -- Extract each frame in turn

[nul, cvnRefs, vnDataSize] = GetFullFileRefs(oStack, sSubs.subs);

vnFrames = cvnRefs{1};

if (numel(vnDataSize) == 3)
   tfData = zeros(vnDataSize([1 3]));

elseif (numel(vnDataSize) == 4)
   tfData = zeros(vnDataSize([1 2 4]));
else
   % - Shouldn't ever happen
   error('FocusStack:UnexpectedError', '*** FocusStack/ExtractSummedFrames: Unexpected error in data sizing.');
end

% - Can we return aligned frames, if requested?
wOld = warning;
if (bAligned && isempty(oStack.mfFrameShifts))
   % - No, so just return un-aligned data
   warning('FocusStack:UnalignedStack', '--- FocusStack/ExtractSummedFrames: Warning: The stack has not yet been aligned.  Returning un-aligned stack images.');
   bAligned = false;
   warning('off', 'FocusStack:UnalignedStack');
end

for (nFrame = vnFrames(:)')
   % - Extract a single frame and accumulate
   if (~bAligned)
      tfThisFrame = oStack.ExtractFrames({cvnRefs{3}, nFrame, cvnRefs{2}});
   else
      tfThisFrame = oStack.ExtractAlignedFrames({cvnRefs{3}, nFrame, cvnRefs{2}});
   end
   tfData = tfData + double(reshape(tfThisFrame, size(tfData)));
end

% - Restore warning state;
warning(wOld);

% --- END of ExtractSummedFrames METHOD ---
