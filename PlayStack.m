function [hFigure] = PlayStack(oStack, vnChannels, vfDataRange)

% PlayStack - FUNCTION Make a window to play or scrub through a stack
%
% Usage: [hFigure] = PlayStack(oStack, vnChannels, vfDataRange)
%        [hFigure] = PlayStack(oStack, fhExtractionFunction, vfDataRange)
%
% 'oStack' is a 4D tensor (X Y nFrame nChannel) or a FocusStack object.
%
% Uses 'videofig' from Joo Filipe Henriques.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>

% -- Check arguments

if (nargin < 1)
   disp('*** PlayStack: Incorrect usage.');
   help PlayStack;
   return;
end

if (~exist('vnChannels', 'var') || isempty(vnChannels))
   vnChannels = 1;
end

if (isa(vnChannels, 'function_handle'))
   if (~isa(oStack, 'FocusStack'))
      disp('*** PlayStack: If an extraction function is provided, a FocusStack object must also be provided.')
      return;
   end
   
   fhExtractionFunc = vnChannels;
else
   if (isa(oStack, 'FocusStack'))
      fhExtractionFunc = @(fsData, vnPixels, vnFrames)fsData.AlignedStack(vnPixels, vnFrames, vnChannels);
   else
      fhExtractionFunc = @(fsData, vnPixels, vnFrames)reshape(fsData(:, :, vnFrames, vnChannels), 1, [], numel(vnChannels));
   end
end

if (~exist('vfDataRange', 'var'))
   vfDataRange = [];
end

% -- Get stack parameters
vnStackSize = size(oStack);
nStackLength = vnStackSize(3);
vnFrameSize = vnStackSize(1:2);

if (isfield(oStack, 'tFrameDuration'))
   tFPS = 1 ./ oStack.tFrameDuration;
else
   tFPS = [];
end

if (isa(oStack, 'FocusStack'))
   mbAlignedMask = oStack.GetAlignedMask();
else
   mbAlignedMask = true(size(oStack, 1), size(oStack, 2));
end

% - Make a videofig and show the first frame
fhRedraw = @(n)(PlotFrame(oStack, n, fhExtractionFunc, vnFrameSize, nStackLength, mbAlignedMask, vfDataRange));
hFigure = videofig(  nStackLength, ...
                     fhRedraw, ...
                     tFPS);
fhRedraw(1);

% - Remove return argument, if not requested
if (nargout == 0)
   clear hFigure;
end


% --- END of PlayStack FUNCTION ---

function PlotFrame(oStack, nFrame, fhExtractionFunc, vnFrameSize, nStackLength, mbDataMask, vfDataRange)

% - Turn off "unaligned stack" warning
wOld = warning('off', 'FocusStack:UnalignedStack');

% - Extract each requested channel
imRGB = zeros([vnFrameSize([2 1]) 3]);

% - Extract frame from the stack
tfFrame = reshape(fhExtractionFunc(oStack, ':', nFrame), vnFrameSize(1), vnFrameSize(2), []);

for (nChannel = 1:size(tfFrame, 3))
   mfThisFrame = tfFrame(:, :, 1, nChannel)';
   
   if (isa(oStack, 'FocusStack') && oStack.bConvertToDFF)
      % - Clip 1..2
      mfThisFrame(mfThisFrame > 2) = 2;
   end
   
   % - Normalise frame within mask
   if (isempty(vfDataRange))
      mfThisFrame = double(mfThisFrame) - min(double(mfThisFrame(mbDataMask)));
      mfThisFrame = uint8(mfThisFrame ./ max(mfThisFrame(mbDataMask)) * 255);
      mfThisFrame(~mbDataMask) = 0;
   else
      mfThisFrame = double(mfThisFrame) - min(vfDataRange);
      mfThisFrame = mfThisFrame ./ abs(diff(vfDataRange));
      mfThisFrame(~mbDataMask) = 0;
      mfThisFrame(mfThisFrame < 0) = 0;
      mfThisFrame(mfThisFrame > 1) = 1;
      mfThisFrame = uint8(mfThisFrame * 255);
   end
   
   % - Assign colour map
   if (nChannel == 1)
      imThisRGB = ind2rgb(mfThisFrame, green);
   elseif (nChannel == 2)
      imThisRGB = ind2rgb(mfThisFrame, red);
   else
      imThisRGB = ind2rgb(mfThisFrame, gray(256));
      imThisRGB = imThisRGB ./ numel(vnChannels);
   end
   
   % - Mix channels
   imRGB = imRGB + imThisRGB;
end

% - Draw image
image(imRGB);
axis tight equal off;

% - Add some information
strTitle = sprintf('Frame [%d] of [%d]', nFrame, nStackLength);
text(5, 5, strTitle, 'Color', 'white');

% - Add a scale bar
if (isa(oStack, 'FocusStack') && ~isempty(oStack.fPixelsPerUM))
   PlotScaleBar(oStack.fPixelsPerUM, 50, 'bl', 'w', 'LineWidth', 6);
end

% - Restore warnings
warning(wOld);


function [mnRedCMap] = red

mnRedCMap = zeros(256, 3);
mnRedCMap(:, 1) = linspace(0, 1, 256);


function [mnGreenCMap] = green

mnGreenCMap = zeros(256, 3);
mnGreenCMap(:, 2) = linspace(0, 1, 256);

% --- END of PlayStack.m ---
