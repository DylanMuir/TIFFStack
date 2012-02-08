function [vfBlackTrace] = DefineBlackRegion(oStack, vnBlackPixels, nChannel)

% DefineBlackRegion - METHOD Define a region of the stack to be "black" for each frame
%
% Subsequent frame extractions will subtract the black value for that frame

if (~exist('nChannel', 'var') || isempty(nChannel))
   nChannel = 1;
end

% - Turn off "unaligned stack" warning
wOld = warning('off', 'FocusStack:UnalignedStack');

% - Disable black subtraction and dF/F conversion
oStack.bSubtractBlack = false;
bConvertToDFF = oStack.bConvertToDFF;
oStack.bConvertToDFF = false;

if (~exist('vnBlackPixels', 'var') || isempty(vnBlackPixels))
   % - Display the first stack frame
   hFigure = figure;
   imagesc(ExtractAlignedFrames(oStack, {':', ':', 1, nChannel})');
   colormap gray;
   axis equal tight off;
   
   % - Let the user define a circular ROI to specify black pixels
   [nul, vnBlackPixels] = circleroi();
   
   % - Transpose X and Y pixels
   vnStackSize = size(oStack);
   [vnY, vnX] = ind2sub(vnStackSize(1:2), vnBlackPixels);
   vnBlackPixels = sub2ind(vnStackSize(1:2), vnX, vnY);
   
   % - Delete the figure
   delete(hFigure);
   drawnow;   
end

% - Extract aligned black trace for the specified pixels
disp('--- FocusStack/DefineBlackRegion: Extracting black trace for specified pixels...');
mfBlackTrace = double(ExtractAlignedFrames(oStack, {vnBlackPixels, ':', nChannel}));

% - Average black trace for each frame
vfBlackTrace = mean(mfBlackTrace, 1);

% - Assign black trace (which turns on black subtraction)
oStack.vfBlackTrace = reshape(vfBlackTrace, [], 1);

% - Restore warnings
warning(wOld);

% - Restore dF/F conversion
oStack.bConvertToDFF = bConvertToDFF;


% - Remove output argument, if not requested
if (nargout == 0)
   clear vfBlackTrace;
end

% --- END of DefineBlackRegion.m ---
