function [hFigure] = PlayStack(oStack, vnChannels)

% PlayStack - FUNCTION Make a window to play or scrub through a stack
%
% Usage: [hFigure] = PlayStack(oStack, vnChannels)
%
% 'oStack' is a 4D tensor (X Y nFrame nChannel)
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

% -- Get stack parameters
vnStackSize = size(oStack);
nStackLength = vnStackSize(3);
vnFrameSize = vnStackSize(1:2);

if (isfield(oStack, 'tFrameDuration'))
   tFPS = 1 ./ oStack.tFrameDuration;
else
   tFPS = [];
end

% - Make a videofig and show the first frame
fhRedraw = @(n)(PlotFrame(oStack, n, vnChannels, vnFrameSize, nStackLength));
hFigure = videofig(  nStackLength, ...
                     fhRedraw, ...
                     tFPS);
fhRedraw(1);

% - Remove return argument, if not requested
if (nargout == 0)
   clear hFigure;
end


% --- END of PlayStack FUNCTION ---

function PlotFrame(oStack, nFrame, vnChannels, vnFrameSize, nStackLength)

% - Turn off "unaligned stack" warning
wOld = warning('off', 'FocusStack:UnalignedStack');

% - Extract each requested channel
imRGB = zeros([vnFrameSize([2 1]) 3]);
for (nChannel = vnChannels(:)')
   % - Extract frame from the stack
   if (isa(oStack, 'FocusStack'))
      mfThisFrame = oStack.AlignedStack(:, :, nFrame, nChannel)';
%       mfThisFrame = oStack.RawStack(:, :, nFrame, nChannel)';
      mbDataMask = oStack.GetAlignedMask';
   else
      mfThisFrame = oStack(:, :, nFrame, nChannel)';
      mbDataMask = true(size(mfThisFrame));
   end
   
   if (isa(oStack, 'FocusStack') && oStack.bConvertToDFF)
      % - Clip 1..2
      mfThisFrame(mfThisFrame > 2) = 2;
   end
   
   % - Normalise frame within mask
   if (isa(mfThisFrame, 'uint8'))
      % - Dont normalise
      
%    elseif (isa(mfThisFrame, 'uint16'))
%       % - Scale to [0..255]
%       mfThisFrame = uint8(double(mfThisFrame) / 256);
      
   else
      mfThisFrame = double(mfThisFrame) - min(double(mfThisFrame(mbDataMask)));
      mfThisFrame = uint8(mfThisFrame ./ max(mfThisFrame(mbDataMask)) * 255);
      mfThisFrame(~mbDataMask) = 0;
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
strTitle = ['Channel(s) [' sprintf('%d,', vnChannels) '], '];
strTitle = [strTitle sprintf('frame [%d] of [%d]', nFrame, nStackLength)];
text(10, 10, strTitle, 'Color', 'white');

% - Add a scale bar
if (isa(oStack, 'FocusStack') && ~isempty(oStack.fPixelsPerUM))
   PlotScaleBar(oStack.fPixelsPerUM, 20, 'bl', 'w', 'LineWidth', 6);
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
