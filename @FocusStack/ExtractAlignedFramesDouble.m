function [tfData] = ExtractAlignedFramesDouble(oStack, sSubs)

% ExtractAlignedFrames - METHOD Extract and align frames from the stack

% - Convert cell references to a subsref structure
if (iscell(sSubs))
   subs = sSubs;
   clear sSubs;
   sSubs.subs = subs;
   sSubs.type = '()';
end

% -- Check whether the reference is even possible

if (numel(sSubs) > 1)
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractFrames: Only one referencing level is possible.');
end

if (~isequal(sSubs.type, '()'))
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractFrames: Only ''()'' referencing is supported.');
end


% -- Check if alignment data exists

if (isempty(oStack.mfFrameShifts))
   % - No, so just return un-aligned data
   warning('FocusStack:UnalignedStack', '--- FocusStack/ExtractAlignedFrames: Warning: The stack has not yet been aligned.  Returning un-aligned stack images.');
   tfData = ExtractFrames(oStack, sSubs);
   return;
end


% -- Create aligned frames cache, if it doesn't exist

vnStackSize = size(oStack);
if (isempty(oStack.vbCachedAlignedFrames) || isempty(oStack.oAlignedFrameCache))
   % - Allocate a cache
   oStack.vbCachedAlignedFrames = false(vnStackSize(3:4));
   oStack.oAlignedFrameCache = MappedTensor(vnStackSize, 'Class', oStack.strDataClass);
end


% -- Align each frame and extract desired pixels

[nul, cvnRefs, vnDataSize] = GetFullFileRefs(oStack, sSubs.subs);

vnFrames = cvnRefs{1};

% - Record data manipulation state
bSubtractBlack = oStack.bSubtractBlack;
bSubtractBlank = oStack.bSubtractBlank;
bConvertToDFF = oStack.bConvertToDFF;

% Always return double (to allow NaN elements for unavailable pixels)
strOutputDataClass = 'double';
% if (bSubtractBlack || bConvertToDFF)
%    strOutputDataClass = 'double';
% else
%    strOutputDataClass = oStack.strDataClass;
% end

% - Determine desired frame size
if (numel(vnDataSize) == 3)
   tfData = zeros(vnDataSize, strOutputDataClass);
   vnFrameSize = vnDataSize(1);

elseif (numel(vnDataSize) == 4)
   tfData = zeros([prod(vnDataSize(1:2)) vnDataSize(3:4)], strOutputDataClass);
   vnFrameSize = vnDataSize(1:2);
   
else
   % - Shouldn't ever happen
   error('FocusStack:UnexpectedError', '*** FocusStack/ExtractAlignedFrames: Unexpected error in data sizing.');
end   

% - Turn off data manipulation
oStack.bSubtractBlack = false;
oStack.bConvertToDFF = false;
oStack.bSubtractBlank = false;
oStack.bDoubleOutput=true;

try
   for (nFrameIndex = 1:numel(vnFrames))
      % -- Check the aligned frame cache
      if (false) % all(oStack.vbCachedAlignedFrames(vnFrames(nFrameIndex), cvnRefs{2})))
         % - Extract the aligned frames from the cache
         tfAlignedFrame = oStack.oAlignedFrameCache(:, :, vnFrames(nFrameIndex), cvnRefs{2});
         
      else
         % -- Extract a full frame from the stack
         % - Permute to get channels in the third dimension
         tfThisFrame = permute(oStack.ExtractFrames({':', ':', vnFrames(nFrameIndex), cvnRefs{2}}), [1 2 4 3]);
         
         % - Compute a translation matrix
         vfThisShift = oStack.mfFrameShifts(vnFrames(nFrameIndex), :);
         
         if (any(vfThisShift ~= 0))
            mfTranslation = eye(3);
            mfTranslation(3, 1:2) = vfThisShift;
            oTForm = maketform('projective', mfTranslation);
            
            % - Align the frame according to the pre-calculated shift data
            tfAlignedFrame = imtransform(tfThisFrame, oTForm, 'XData', [1 vnStackSize(1)], 'YData', [1 vnStackSize(2)], 'FillValues', nan);
         else
            tfAlignedFrame = tfThisFrame;
         end
         
         % - Insert aligned frame into stack
         oStack.oAlignedFrameCache(:, :, vnFrames(nFrameIndex), cvnRefs{2}) = tfAlignedFrame;
         oStack.vbCachedAlignedFrames(vnFrames(nFrameIndex), cvnRefs{2}) = true;
      end
      
      % - Return only the requested pixels
      tfAlignedFrame = reshape(tfAlignedFrame, [], 1, numel(cvnRefs{2}));
      
      % - Should we subtract the black?
      if (bSubtractBlack)
         % - Subtract the black trace for this frame
         tfAlignedFrame = double(tfAlignedFrame) - oStack.vfBlackTrace(vnFrames(nFrameIndex));
      end
      
      % - Reshape data and return
      tfData(:, nFrameIndex, :) = reshape(ipermute(tfAlignedFrame(cvnRefs{3}, :, :), [1 2 4 3]), [prod(vnFrameSize) 1 numel(cvnRefs{2})]);
   end

   % - Should we convert to DFF?
   if (bConvertToDFF)
      % - Extract the blank frame for these frames and divide
      tfBlanks = ExtractBlankFrames(oStack, {cvnRefs{3} cvnRefs{1}});
      tfData = double(tfData) ./ tfBlanks - 1;
   end
   
   % - Should we subtract the blank?
   if (bSubtractBlank)
      % - Extract the blank frame for these frames and subtract
      tfBlanks = ExtractBlankFrames(oStack, {cvnRefs{3} cvnRefs{1}});
      tfData = double(tfData) - tfBlanks;
   end
   
catch err
   % - Restore data manipulation state
   oStack.bSubtractBlack = bSubtractBlack;
   oStack.bConvertToDFF = bConvertToDFF;
   oStack.bSubtractBlank = bSubtractBlank;
   rethrow(err);
end

% - Restore data manipulation state
oStack.bSubtractBlack = bSubtractBlack;
oStack.bSubtractBlank = bSubtractBlank;
oStack.bConvertToDFF = bConvertToDFF;
oStack.bDoubleOutput=false;

% - Return pixels to correct size
tfData = reshape(tfData, vnDataSize);

% --- END ExtractAlignedFrames METHOD ---
