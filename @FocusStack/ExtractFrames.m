function tfData = ExtractFrames(oStack, sSubs)

% ExtractFrames - METHOD Extract frames from files on disk

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


% -- Check referencing and convert linear frame references to file/frame
% references

[cvnFileRefs, cvnFullRefs, vnDataSize] = GetFullFileRefs(oStack, sSubs.subs);


% -- Read from each file in turn

if (oStack.bSubtractBlack || oStack.bConvertToDFF || oStack.bSubtractBlank || oStack.bDoubleOutput)
   strOutputDataClass = 'double';
else
   strOutputDataClass = oStack.strDataClass;
end   

% - Pre-allocate data matrix
%   'tfData' is [vnFrames vnChannels vnLinearPixels]
tfData = zeros(numel(cvnFileRefs{2}), numel(cvnFileRefs{3}), numel(cvnFileRefs{4}), strOutputDataClass);

vnUniqueFiles = unique(cvnFileRefs{1});
for (nFile = vnUniqueFiles(:)')
   vbMatchingRefs = cvnFileRefs{1} == nFile;
   tfPixels = ReadPixelsFromFile(oStack, nFile, cvnFileRefs{4}, cvnFileRefs{2}(vbMatchingRefs), cvnFileRefs{3});
   
   % - Should we subtract the black?
   if (oStack.bSubtractBlack)
      % - Extract the black trace for these pixels
      vfThisBlackTrace = oStack.vfBlackTrace(cvnFullRefs{1}(vbMatchingRefs));
      tfPixels = double(tfPixels) - repmat(vfThisBlackTrace, [1 numel(cvnFullRefs{2}) numel(cvnFullRefs{3})]);
   end

   % - Should we convert to DFF?
   if (oStack.bConvertToDFF)
      % - Extract the blank pixels for these pixels and divide
      tfBlankPixels = permute(ExtractBlankFrames(oStack, [cvnFullRefs(3) {cvnFullRefs{1}(vbMatchingRefs)}]), [2 3 1]);
      tfPixels = double(tfPixels) ./ tfBlankPixels - 1;
   end
   
   % - Should we subtract the blank?
   if (oStack.bSubtractBlank)
      tfBlankPixels = permute(ExtractBlankFrames(oStack, [cvnFullRefs(3) {cvnFullRefs{1}(vbMatchingRefs)}]), [2 3 1]);
      tfPixels = double(tfPixels) - tfBlankPixels;
   end
   
   % - Assign data to be returned
   tfData(vbMatchingRefs, :, :) = tfPixels;
end

% - Re-shape data
tfData = reshape(permute(tfData, [3 1 2]), vnDataSize);


function tfData = ReadPixelsFromFile(oStack, nFileNumber, vnPixels, vnFrames, vnChannels)

% - tfData is [vnFrames vnChannels vnPixels]

switch (class(oStack.vhMemMapFileHandles{nFileNumber}))
   case 'memmapfile'
      try
         tfData = permute(oStack.vhMemMapFileHandles{nFileNumber}.Data.tfStack(vnChannels, vnPixels, vnFrames), ...
                          [3 1 2]);
      catch err
         disp('meh');
         rethrow(err);
      end

   case 'MappedTensor'
      tfData = permute(oStack.vhMemMapFileHandles{nFileNumber}(vnChannels, vnPixels, vnFrames), [3 1 2]);
      
      % - Get full frames and subsample
%       tfData = permute(oStack.vhMemMapFileHandles{nFileNumber}(vnChannels, :, vnFrames), [3 1 2]);
%       tfData = tfData(:, :, vnPixels);

   case 'TIFFStack'
      % - Get full frames
      tfData = oStack.vhMemMapFileHandles{nFileNumber}(:, :, vnFrames, vnChannels);
      
      % - Subsample frames
      tfData = reshape(tfData, 1, [], numel(vnFrames), numel(vnChannels));
      tfData = permute(tfData(1, vnPixels, :, :), [3 4 2 1]);
      
   otherwise
      error('FocusStack:UnexpectedError', '*** FocusStack/ExtractFrames: Unexpected stack type');
end


% --- END of ExtractFrames METHOD ---
