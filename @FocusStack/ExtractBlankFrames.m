function [tfBlankMean, tfBlankStd] = ExtractBlankFrames(oStack, sSubs)

% ExtractBlankFrames - METHOD Extract blank frames (pixels) from a stack
%
% Usage: [tfBlankMean, tfBlankStd] = ExtractBlankFrames(oStack, sSubs)
%
% ExtractBlankFrames just pulls pixel traces from the stack corresponding to the
% referenced pixels, but with the pixel data replaced with the blank frame data
% (and blank standard deviation data).
%
% This can use a lot of memory.  You might find the method
% GetCorrespondingBlankFrames more useful.

% - Convert cell references to a subsref structure
if (iscell(sSubs))
   subs = sSubs;
   clear sSubs;
   sSubs.subs = subs;
   sSubs.type = '()';
end


% -- Check whether the reference is even possible

if (numel(sSubs) > 1)
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractBlankFrames: Only one referencing level is possible.');
end

if (numel(sSubs.subs) > 3)
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractBlankFrames: No channel reference is allowed for blank frames.');
end

if (~isequal(sSubs.type, '()'))
   error('FocusStack:ImproperReference', '*** FocusStack/ExtractBlankFrames: Only ''()'' referencing is supported.');
end


% -- Check referencing and convert linear frame references to file/frame
% references

% - Pretend we're looking at channel '1'
cSubs = sSubs.subs;
cSubs{end+1} = 1;

[nul, cvnFullRefs, vnDataSize] = GetFullFileRefs(oStack, cSubs);

bGetStds = nargout > 1;

if (prod(vnDataSize) == 0)
    tfBlankMean = nan(vnDataSize);

    if (bGetStds)
        tfBlankStd = nan(vnDataSize);
    end
    
    return;
end

% -- Extract unique blank frames

[tfUniqueBlankMeans, tfUniqueBlankStds, vnBlankFrameIndices] = ...
   GetCorrespondingBlankFrames(oStack, sSubs);

nNumUniqueBlanks = size(tfUniqueBlankMeans, 3);


% -- Replicate blank pixels for all frames

tfBlankMean = ones([numel(cvnFullRefs{3}) numel(cvnFullRefs{1}) numel(cvnFullRefs{2})]);

if (bGetStds)
   tfBlankStd = ones([numel(cvnFullRefs{3}) numel(cvnFullRefs{1})]);
end

for (nUniqueBlankIndex = 1:nNumUniqueBlanks)
   % - Find stack frames that use this blank frame
   vbMatchingFrames = vnBlankFrameIndices == nUniqueBlankIndex;
   
   % - Assign blank pixels
   tfBlankMean(:, vbMatchingFrames) = repmat(reshape(tfUniqueBlankMeans(:, :, nUniqueBlankIndex), [], 1), 1, nnz(vbMatchingFrames));
   
   if (bGetStds)
      tfBlankStd(:, vbMatchingFrames) = repmat(reshape(tfUniqueBlankStds(:, :, nUniqueBlankIndex), [], 1), 1, nnz(vbMatchingFrames));
   end
end

% - Reshape blank frames
tfBlankMean = reshape(tfBlankMean, vnDataSize);

if (bGetStds)
   tfBlankStd = reshape(tfBlankStd, vnDataSize);
end

% --- END of ExtractBlankFrames METHOD ---
