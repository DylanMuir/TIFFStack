function [tfBlankFrames, tfBlankStds, vnCorrespondingBlanks] = ...
   GetCorrespondingBlankFrames(oStack, sSubs)

% GetCorrespondingBlankFrames - METHOD Get the blank frames corresponding to a set of pixels
%
% Usage: [tfBlankFrames, tfBlankStds, vnCorrespondingBlanks] = ...
%    GetCorrespondingBlankFrames(oStack, sSubs)
%
% 'oStack' is a FocusStack.  'sSubs' is a subsref structure referencing the
% stack.
%
% 'tfBlankFrames' and 'tfBlankStds' will be a set of blank frames, [UxVxN],
% where UxV is the dimensions of the pixels accessed from each frame in 'sSubs'.
% This could also be linear indexing, in which case it will be [Ux1xN].  'N' is
% NOT the total number of frames indexed in 'sSubs'.  Instead,
% 'vnCorrespondingBlanks' indicates for each frame in 'sSubs', which slice of
% 'tfBlankFrames' and 'tfBlankStds' corresponds to that frame.

% - Convert cell references to a subsref structure
if (iscell(sSubs))
   subs = sSubs;
   clear sSubs;
   sSubs.subs = subs;
   sSubs.type = '()';
end


% -- Check whether the reference is even possible

if (numel(sSubs) > 1)
   error('FocusStack:ImproperReference', '*** FocusStack/GetCorrespondingBlankFrames: Only one referencing level is possible.');
end

if (numel(sSubs.subs) > 3)
   error('FocusStack:ImproperReference', '*** FocusStack/GetCorrespondingBlankFrames: No channel reference is allowed for blank frames.');
end

if (~isequal(sSubs.type, '()'))
   error('FocusStack:ImproperReference', '*** FocusStack/GetCorrespondingBlankFrames: Only ''()'' referencing is supported.');
end


% -- Check referencing and convert linear frame references to file/frame
% references

% - Pretend we're looking at channel 1
cSubs = sSubs.subs;
cSubs{end+1} = 1;

[nul, cvnFullRefs, vnDataSize] = GetFullFileRefs(oStack, cSubs);

if (numel(vnDataSize) == 3)
   vnOutDataSize = [vnDataSize(1) 1 1 vnDataSize(3)];
   nNumFrames = vnDataSize(2);
else
   vnOutDataSize = [vnDataSize(1:2) 1 vnDataSize(4)];
   nNumFrames = vnDataSize(3);
end


% -- Extract blank frames indices for each frame

bGetMeans = ~isempty(oStack.mnAssignedBlankMeanFrames);
bGetStds = ~isempty(oStack.mnAssignedBlankStdFrames) && (nargout > 1);

% - Have any blank frames been assigned?
if (~bGetMeans)
   % - No, so just return "one" pixels
   tfBlankFrames = ones(vnOutDataSize);
   vnCorrespondingBlanks = ones(nNumFrames, 1);
   
   if (nargout > 1)
      tfBlankStds = ones(vnOutDataSize);
   end
   
   return;
end

[mnFrames, mnChannels] = ndgrid(cvnFullRefs{1}, cvnFullRefs{2});

if (bGetMeans)
   vnBlankFrameMeanCacheIndices = sub2ind(size(oStack.mnAssignedBlankMeanFrames), mnFrames(:), mnChannels(:));
   vnBlankFrameMeanIndices = oStack.mnAssignedBlankMeanFrames(vnBlankFrameMeanCacheIndices);
end

if (bGetStds)
   vnBlankFrameStdCacheIndices = sub2ind(size(oStack.mnAssignedBlankStdFrames), mnFrames(:), mnChannels(:));
   vnBlankFrameStdIndices = oStack.mnAssignedBlankStdFrames(vnBlankFrameStdCacheIndices);
end

% - Check to see if blank frames exist
if ((bGetMeans && any(isnan(vnBlankFrameMeanIndices))) || ...
      bGetStds && any(isnan(vnBlankFrameStdIndices)))
   error('FocusStack:UnassignedBlank', ...
      '*** FocusStack/GetCorrespondingBlankFrames: Not all requested stack frames have an associated blank.');
end


% -- Extract required blank frames and pull out required pixels

if (bGetMeans)
   [vnUniqueBlankMeanFrames, nul, vnCorrespondingBlanks] = unique(vnBlankFrameMeanIndices);
   cmfTheseMeanBlanks = oStack.cmfBlankFrames(vnUniqueBlankMeanFrames);

%    sSubs = substruct('()', cvnFullRefs(3));
%    cmfTheseMeanBlanks = cellfun(@(c)subsref(c, sSubs), cmfTheseMeanBlanks, 'UniformOutput', false);
%    cmfTheseMeanBlanks = ref_cells(cmfTheseMeanBlanks, cvnFullRefs{3});
   tfBlankFrames = reshape(cat(3, cmfTheseMeanBlanks{:}), [], 1, numel(vnUniqueBlankMeanFrames));
   tfBlankFrames = tfBlankFrames(cvnFullRefs{3}, :, :);
   tfBlankFrames = reshape(tfBlankFrames, vnOutDataSize(1), vnOutDataSize(2), [], vnOutDataSize(4));
end

if (bGetStds)
   vnUniqueBlankStdFrames = unique(vnBlankFrameStdIndices);
   cmfTheseStdBlanks = oStack.cmfBlankFrames(vnUniqueBlankStdFrames);

%    sSubs = substruct('()', cvnFullRefs(3));
%    cmfTheseStdBlanks = cellfun(@(c)subsref(c, sSubs), cmfTheseStdBlanks, 'UniformOutput', false);
%    cmfTheseStdBlanks = ref_cells(cmfTheseStdBlanks, cvnFullRefs{3});
%    tfBlankStds = reshape(cat(3, cmfTheseStdBlanks{:}), vnOutDataSize(1), vnOutDataSize(2), [], vnOutDataSize(4));
   tfBlankStds = reshape(cat(3, cmfTheseStdBlanks{:}), [], 1, numel(vnUniqueBlankMeanFrames));
   tfBlankStds = tfBlankStds(cvnFullRefs{3}, :, :);
   tfBlankStds = reshape(tfBlankStds, vnOutDataSize(1), vnOutDataSize(2), [], vnOutDataSize(4));
end


function cellFrames = ref_cells(cellFrames, vnPixels) %#ok<DEFNU>

for (nCellIndex = 1:numel(cellFrames)) %#ok<FORPF>
   cellFrames{nCellIndex} = cellFrames{nCellIndex}(vnPixels);
end

% --- END of GetCorrespondingBlankFrames.m ---
