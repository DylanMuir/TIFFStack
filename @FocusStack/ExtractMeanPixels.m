function [tfData] = ExtractMeanPixels(oStack, sSubs, bAligned)

% ExtractMeanPixels - METHOD Extract averaged pixels from a stack

% -- Check whether the reference is even possible

if (numel(sSubs) > 1)
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractMeanPixels: Only one referencing level is possible.');
end

if (~isequal(sSubs.type, '()'))
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractMeanPixels: Only ''()'' referencing is supported.');
end


% -- Check arguments

if (~exist('bAligned', 'var') || isempty(bAligned))
   bAligned = false;
end


% -- Extract each frame in turn

[nul, cvnRefs, vnDataSize] = GetFullFileRefs(oStack, sSubs.subs);

vnFrames = cvnRefs{1};

if (numel(vnDataSize) == 3)
   tfData = zeros([1 vnDataSize([2 3])]);

elseif (numel(vnDataSize) == 4)
   tfData = zeros([1 vnDataSize([3 4])]);
   
else
   % - Shouldn't ever happen
   error('FocusStack:UnexpectedError', '*** FocusStack/ExtractMeanPixels: Unexpected error in data sizing.');
end

% - Can we return aligned frames, if requested?
wOld = warning;
if (bAligned && isempty(oStack.mfFrameShifts))
   % - No, so just return un-aligned data
   warning('FocusStack:UnalignedStack', '--- FocusStack/ExtractMeanPixels: Warning: The stack has not yet been aligned.  Returning un-aligned stack images.');
   bAligned = false;
   warning('off', 'FocusStack:UnalignedStack');
end

for (nFrame = 1:numel(vnFrames))
   % - Extract a single frame and accumulate
   if (~bAligned)
      tfThisFrame = oStack.ExtractFrames({cvnRefs{3}, vnFrames(nFrame), cvnRefs{2}});
   else
      tfThisFrame = oStack.ExtractAlignedFrames({cvnRefs{3}, vnFrames(nFrame), cvnRefs{2}});
   end
   
   tfData(1, nFrame, :) = mean(tfThisFrame, 1);
end

% - Restore warning state;
warning(wOld);

% --- END of ExtractMeanPixels METHOD ---
