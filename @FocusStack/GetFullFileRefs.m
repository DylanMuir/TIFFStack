function [cvnFileRefs, cvnFullRefs, vnDataSize] = GetFullFileRefs(oStack, cRefs)

% GetFullFileRefs - PRIVATE METHOD Check and convert references
%
% Usage: [cvnFileRefs, cvnFullRefs, vnDataSize] = GetFullFileRefs(oStack, cRefs)
%
% 'cvnFileRefs' is [vnFiles vnFramesInFile vnChannels vnPixels]
% 'cvnFullRefs' is [vnFrames vnChannels vnPixels]
% 'vnDataSize' is the desired return size of the extracted pixels

vnStackSize = size(oStack);
nNumRefDims = numel(cRefs);

% -- Deal with linear or matrix within-frame references

% if (nNumRefDims == 1)
%    % -- Full linear indexing
%    % - Return linear data
%    vnStackSize = [prod(vnStackSize) 1];
%    cRefs = ConvertColonsCheckLims(cRefs, vnStackSize);
%    
%    % - Convert to frame-linear/frame/channel refs
%    [cvnRefs{1:3}] = ind2sub([prod(vnStackSize(1:2)) vnStackSize(3:4)], cRefs{:});
%    vnDataSize = size(cRefs{1});
   
if (nNumRefDims == 3)
   % -- Frame-linear indexing
   % - Return frame-linear/frame/channel data
   vnStackSize = [prod(vnStackSize(1:2)) vnStackSize(3:4)];
   cvnRefs = ConvertColonsCheckLims(cRefs, vnStackSize);
   
   % - Use frame-linear/frame/channel refs
   vnDataSize = [numel(cvnRefs{1}) numel(cvnRefs{2}) numel(cvnRefs{3})];
   
elseif (nNumRefDims == 4)
   % -- Matrix-frame indexing
   % - Return frame-matrix/frame/channel data
   cRefs = ConvertColonsCheckLims(cRefs, vnStackSize);
   
   % - Convert to frame-linear/frame/channel refs
   [vnRows, vnCols] = ndgrid(cRefs{1}, cRefs{2});
   [cvnRefs{1}] = sub2ind(vnStackSize(1:2), vnRows(:), vnCols(:));
   cvnRefs(2:3) = cRefs(3:4);
   vnDataSize = [numel(cRefs{1}) numel(cRefs{2}) numel(cRefs{3}) numel(cRefs{4})];

else
   % -- Invalid reference
   error('FocusStack:InvalidRef', '*** FocusStack/GetFullFileRefs: Invalid reference format.');
end



% -- Determine which file(s) and which frames are needed

vnFrames = reshape(cvnRefs{2}, 1, []);
vnFiles = ones(1, numel(vnFrames));

bAdjustFiles = true;
while (bAdjustFiles)
   vbOverflow = vnFrames > oStack.vnNumFrames(vnFiles);
   vnFrames(vbOverflow) = vnFrames(vbOverflow) - oStack.vnNumFrames(vnFiles(vbOverflow));
   vnFiles(vbOverflow) = vnFiles(vbOverflow) + 1;

   % - Do we need to continue?
   bAdjustFiles = any(vbOverflow);
end

cvnFileRefs = {vnFiles vnFrames cvnRefs{3} cvnRefs{1}};
cvnFullRefs = cvnRefs([2 3 1]);


function [cCheckedRefs] = ConvertColonsCheckLims(cRefs, vnLims)

% - Check each dimension in turn
for (nRefDim = numel(cRefs):-1:1) %#ok<FORPF>
   % - Convert colon references
   if (ischar(cRefs{nRefDim}) && isequal(cRefs{nRefDim}, ':'))
      cCheckedRefs{nRefDim} = 1:vnLims(nRefDim); %#ok<AGROW>
   
   elseif (islogical(cRefs{nRefDim}))
      % - Logical referencing -- convert to linear referencing
      vnIndices = find(cRefs{nRefDim}(:));
      if (any(vnIndices > vnLims(nRefDim)))
         error('FocusStack:InvalidRef', ...
            '*** FocusStack/GetFullFileRefs: Logical referencing for dimension [%d] was out of bounds [1..%d].', ...
            nRefDim, vnLims(nRefDim));
      end      
      cCheckedRefs{nRefDim} = vnIndices; %#ok<AGROW>
      
   elseif (any(cRefs{nRefDim}(:) < 1) || any(cRefs{nRefDim}(:) > vnLims(nRefDim)))
      % - Check limits   
      error('FocusStack:InvalidRef', ...
         '*** FocusStack/GetFullFileRefs: Reference dimension [%d] was out of bounds [1..%d].', ...
         nRefDim, vnLims(nRefDim));
   
   else
      % - This dimension was ok
      cCheckedRefs{nRefDim} = cRefs{nRefDim}(:); %#ok<AGROW>
   end
end


% --- END of GetFullFileRefs METHOD ---

