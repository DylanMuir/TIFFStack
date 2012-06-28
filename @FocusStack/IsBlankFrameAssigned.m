function [vbBlankMeanAssigned, vbBlankStdAssigned] = IsBlankFrameAssigned(oStack, vnFrames)

% IsBlankFrameAssigned - METHOD Queries whether a blank frame has been assigned for the corresponding frames
%
% Usage: [vbBlankMeanAssigned, vbBlankStdAssigned] = IsBlankFrameAssigned(oStack, vnFrames)
%
% 'oStack' is a FocusStack.  'vnFrames' is a vector of frame indices.
%
% 'vbBlankAssigned' will be a boolean vector the same size as 'vnFrames',
% indicating for each corresponding entry in 'vnFrames' whether a blank
% frame was assigned for that frame.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 28th June, 2012

% -- Check whether the reference is even possible

if (any(vnFrames(:) < 1) || any(vnFrames(:) > size(oStack, 3)))
   error('FocusStack:ImproperReference', '*** FocusStack/IsBlankFrameAssigned: Index exceeds matrix dimensions.');
end


vnBlankFrameMeanCacheIndices = sub2ind(size(oStack.mnAssignedBlankMeanFrames), vnFrames(:), 1);
vnBlankFrameMeanIndices = oStack.mnAssignedBlankMeanFrames(vnBlankFrameMeanCacheIndices);

vnBlankFrameStdCacheIndices = sub2ind(size(oStack.mnAssignedBlankStdFrames), vnFrames(:), 1);
vnBlankFrameStdIndices = oStack.mnAssignedBlankStdFrames(vnBlankFrameStdCacheIndices);

vbBlankMeanAssigned = ~isnan(vnBlankFrameMeanIndices);
vbBlankStdAssigned = ~isnan(vnBlankFrameStdIndices);

% --- END of IsBlankFrameAssigned.m ---
