function [mbAlignedMask] = GetAlignedMask(oStack)

% - GetAlignedMask - METHOD Extract a mask for valid aligned data

% - Pre-allocate the mask
vnStackSize = size(oStack);
mbAlignedMask = true(vnStackSize(1:2));

% - Is the stack aligned?
if (isempty(oStack.mfFrameShifts))
   % - No, so just return the "unmasked" mask
   return;
end

% - The stack has been aligned, so work out which pixels are valid

vnMaxNegShift = floor(min(oStack.mfFrameShifts));
vnMaxPosShift = ceil(max(oStack.mfFrameShifts));

mbAlignedMask(1:vnMaxPosShift(2), :) = false;
mbAlignedMask(:, 1:vnMaxPosShift(1)) = false;

mbAlignedMask((end+vnMaxNegShift(2)+1):end, :) = false;
mbAlignedMask(:, (end+vnMaxNegShift(1)+1):end) = false;

% --- END of GetAlignedMask METHOD ---
