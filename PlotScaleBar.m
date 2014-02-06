function [hBar] = PlotScaleBar(fPixPerMM, fLengthMM, strPos, varargin)

% PlotScaleBar - FUNCTION Add a scale bar to the current plot
% 
% Usage: [hBar] = PlotScaleBar(fPixPerMM <, fLengthMM, strPos, PlotOptions...>)

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 30th August, 2006

% -- Default arguments

DEF_fLengthMM = 1;
DEF_strPos = 'tr';
DEF_PlotOptions = {'k-', 'LineWidth', 5};

fRightPos = 0.05;
fLeftPos = 0.95;
fTopPos = 0.05;
fBottomPos = 0.95;


% -- Check arguments

if (nargin < 1)
   disp('*** PlotScaleBar: Incorrect usage')
   help PlotScaleBar;
   return;
end


% -- Set defaults

if (~exist('fLengthMM', 'var') || isempty(fLengthMM))
   fLengthMM = DEF_fLengthMM;
end

if (~exist('strPos', 'var') || isempty(strPos))
   strPos = DEF_strPos;
end

if (nargin > 3)
   PlotOptions = varargin;
else
   PlotOptions = DEF_PlotOptions;
end


% -- Get axis limits

vAxis = axis;

if (strcmp(get(gca, 'XDir'), 'reverse'))
   vAxis([1 2]) = vAxis([2 1]);
end

if (strcmp(get(gca, 'YDir'), 'reverse'))
   vAxis([3 4]) = vAxis([4 3]);
end


% -- How long is the scale bar?

fBarLengthPix = fLengthMM * fPixPerMM;


% -- Where should the scale bar be plotted?

switch lower(strPos)
   case {'tr', 'rt'}
      fWidth = diff(vAxis([1 2]));
      fEndX = vAxis(2) - fWidth * fRightPos;
      fStartX = fEndX - fBarLengthPix * sign(fWidth);
      
      fEndY = vAxis(4) - diff(vAxis([3 4])) * fTopPos;
      fStartY = fEndY;
      fTextY = fStartY + 5;
      
   case {'tl', 'lt'}
      fWidth = diff(vAxis([1 2]));
      fEndX = vAxis(2) - fWidth * fLeftPos;
      fStartX = fEndX + fBarLengthPix * sign(fWidth);
      
      fEndY = vAxis(4) - diff(vAxis([3 4])) * fTopPos;
      fStartY = fEndY;
      fTextY = fStartY + 5;

   case {'br', 'rb'}
      fWidth = diff(vAxis([1 2]));
      fEndX = vAxis(2) - fWidth * fRightPos;
      fStartX = fEndX - fBarLengthPix * sign(fWidth);
      
      fEndY = vAxis(4) - diff(vAxis([3 4])) * fBottomPos;
      fStartY = fEndY;
      fTextY = fStartY - 5;
      
   case {'bl', 'lb'}
      fWidth = diff(vAxis([1 2]));
      fEndX = vAxis(2) - fWidth * fLeftPos;
      fStartX = fEndX + fBarLengthPix * sign(fWidth);
      
      fEndY = vAxis(4) - diff(vAxis([3 4])) * fBottomPos;
      fStartY = fEndY;
      fTextY = fStartY - 5;

   otherwise
      disp('*** PlotScaleBar: Unknown position option.  Must be a combination of l/r and t/b');
      return;
end


% -- Plot line

bIsHold = ishold;
hold on;
hBar = plot([fStartX fEndX], [fStartY fEndY], PlotOptions{:});
text(mean([fStartX fEndX]), fTextY, sprintf('%d um', round(fLengthMM * 1e3)), ...
   'FontSize', 18, 'HorizontalAlignment', 'center', 'Color', get(hBar, 'Color'));

% - Restore hold value
if (~bIsHold)
   hold off;
end

% - Clear return handle, if not requested
if (nargout == 0)
   clear hBar;
end

% --- END of PlotScaleBar.m ---
