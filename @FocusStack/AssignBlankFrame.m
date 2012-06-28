function [mfBlankFrame] = AssignBlankFrame(oStack, mfBlankFrame, vnStackFrames, vnChannels)

% AssignBlankFrame - METHOD Assign a blank frame to a stretch of stack frames, for use in DFF conversion
%
% Usage: AssignBlankFrame(oStack, mfBlankFrame <, vnStackFrames, vnChannels>)
%        AssignBlankFrame(oStack, tfBlankFrame <, vnStackFrames, vnChannels>)
%
% Asign a blank frame (with an optional standard deviation frame) to a set
% of frames in a stack.  If 'tfBlankFrame' is supplied as an [NxMx2]
% tensor, then the second slice is taken as the standard deviation frame.

% -- Check arguments

if (nargin < 2)
   error('FocusStack:BadArguments', '*** FocusStack/AssignBlankFrame: Incorrect usage.');
end

vnStackSize = size(oStack);

vnBlankSize = size(mfBlankFrame);
if (~isequal(vnBlankSize(1:2), vnStackSize(1:2)))
   error('FocusStack:BadBlankFrame', ...
      '*** FocusStack/AssignBlankFrame: ''mfBlankFrame'' must be size [%d %d] for this stack.', vnStackSize(1:2));
end

if (~exist('vnStackFrames', 'var') || isempty(vnStackFrames))
   vnStackFrames = 1:vnStackSize(3);
end

if (~exist('vnChannels', 'var') || isempty(vnChannels))
   vnChannels = 1;
end

% - Convert to numerical indexing
if (islogical(vnStackFrames))
   vnStackFrames = find(vnStackFrames);
end

if (any(vnStackFrames(:) < 1) || any(vnStackFrames(:) > vnStackSize(3)))
   error('FocusStack:BadStackFrameID', ...
      '*** FocusStack/AssignBlankFrame: ''vnStackFrames'' must be limited to [1..%d] for this stack.', vnStackSize(3));
end   

if (any(vnChannels(:) < 1) || any(vnChannels(:) > vnStackSize(4)))
   error('FocusStack:BadChannelID', ...
      '*** FocusStack/AssignBlankFrame: ''vnChannels'' must be limited to [1..%d] for this stack.', vnStackSize(4));
end

if (size(mfBlankFrame, 3) > 1)
   mfStdFrame = mfBlankFrame(:, :, 2);
   mfBlankFrame = mfBlankFrame(:, :, 1);
end


% -- Get a new blank frame ID

nBlankFrameID = numel(oStack.cmfBlankFrames) + 1;

if (isempty(oStack.mnAssignedBlankMeanFrames))
   oStack.mnAssignedBlankMeanFrames = nan(vnStackSize(3:4));
   oStack.mnAssignedBlankStdFrames = nan(vnStackSize(3:4));
end

% - Filter blank frame
if (any(mfBlankFrame(:) < 2))
   warning('FocusStack:FilteringBlankFrames', ...
      '--- FocusStack/AssignBlankFrames: Warning: This blank frame required filtering.');
   
   % - Filter by clipping small values
   mfBlankFrame(mfBlankFrame < 2) = 2;
end

if (any(isnan(mfBlankFrame(:))) || (exist('mfStdFrame', 'var') && any(isnan(mfStdFrame(:)))))
   warning('FocusStack:BlankContainsNaN', ...
      '--- FocusStack/AssignBlankFrames: Warning: This blank frame contains NaNs, which can corrupt further processing.');
end

% - Assign blank mean frame and associate with stack frames
oStack.cmfBlankFrames{nBlankFrameID} = double(mfBlankFrame);
oStack.mnAssignedBlankMeanFrames(vnStackFrames(:), vnChannels(:)) = nBlankFrameID;

% - Filter and assign blank standard deviation frame, if it exists
if (exist('mfStdFrame', 'var'))
   if (any(mfStdFrame(:) == 0))
      warning('FocusStack:FilteringBlankFrames', ...
         '--- FocusStack/AssignBlankFrames: Warning: This blank std. dev. frame required filtering.');
      
      % - Filter by making zero std. devs. equal to the mean std. dev.
      mfStdFrame(mfStdFrame == 0) = mean(mfStdFrame(:));
   end
   
   oStack.cmfBlankFrames{nBlankFrameID+1} = double(mfStdFrame);
   oStack.mnAssignedBlankStdFrames(vnStackFrames(:), vnChannels(:)) = nBlankFrameID+1;
end

% --- END of AssignBlankFrame.m ---
