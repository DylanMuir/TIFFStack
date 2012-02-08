function [hFigure] = MakeROIOverviewFigure(fsData, fPixPerUM, vnChannels, sCellRegions, vnLabelRegions, bSimplePlot)

% MakeROIOverviewFigure - FUNCTION
%
% Usage: [hFigure] = MakeROIOverviewFigure(fsData <, fPixPerUM, vnChannels, sCellRegions, vnLabelRegions, bSimplePlot>)

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 23rd May, 2011

% -- Defaults

DEF_vnChannels = 1:size(fsData, 4);


% -- Check arguments

if (nargin < 1)
   disp('*** MakeROIOverviewFigure: Incorrect usage');
   help MakeROIOverviewFigure;
   return;
end

if (~exist('vnChannels', 'var') || isempty(vnChannels))
   vnChannels = DEF_vnChannels;
end

if (~exist('sCellRegions', 'var'))
   bPlotRegions = false;
else
   bPlotRegions = true;
end

if ((~exist('vnLabelRegions', 'var') || isempty(vnLabelRegions)) && bPlotRegions)
   vnLabelRegions = 1:sCellRegions.NumObjects;
end

if (~exist('bSimplePlot', 'var') || isempty(bSimplePlot))
   bSimplePlot = false;
end

if (isa(fsData, 'double'))
   tfAvgSignal = fsData;
   vnStackSize = size(tfAvgSignal);
   vnStackSize = vnStackSize(1:2);
   fPixPerUM = nan;
else
   % - Get average frames
   strOldNorm = fsData.BlankNormalisation('none');
   tfAvgSignal = fsData.SummedAlignedFrames(:, :, :, vnChannels);
   vnStackSize = size(fsData, 1:2);

   % - Restore normalisation
   fsData.BlankNormalisation(strOldNorm);
end

if (isempty(fPixPerUM))
   fPixPerUM = fsData.fPixelsPerUM;
end


% -- Make a figure

imAvg = tfAvgSignal;

if (numel(vnChannels) == 1)
   tfAvgSignal(:, :, 2) = 0;
end

% - Transpose image and normalise
imAvg(:, :, 1) = tfAvgSignal(:, :, 2)' - min(min(tfAvgSignal(:, :, 2)));
imAvg(:, :, 1) = imAvg(:, :, 1) ./ max(max(imAvg(:, :, 1)));
imAvg(:, :, 2) = tfAvgSignal(:, :, 1)' - min(min(tfAvgSignal(:, :, 1)));
imAvg(:, :, 2) = imAvg(:, :, 2) ./ max(max(imAvg(:, :, 2)));
imAvg(:, :, 3) = 0;

newplot;
hFigure = gcf;

if (~bSimplePlot)
   subplot(1, 3, 1);
   imagesc(tfAvgSignal(:, :, 1)');
   axis equal tight off;
   colormap gray;
   PlotScaleBar(fPixPerUM, 25, 'tr', 'w-', 'LineWidth', 8);
   
   subplot(1, 3, 2);
   image(imAvg);
   axis equal tight off;
   PlotScaleBar(fPixPerUM, 25, 'tr', 'w-', 'LineWidth', 8);
   
   subplot(1, 3, 3);
end

image(imAvg);
axis equal tight off;
hold on;

sContourRegions = sCellRegions;
sContourRegions.NumObjects = nnz(vnLabelRegions);
sContourRegions.PixelIdxList = sContourRegions.PixelIdxList(vnLabelRegions);

contour((labelmatrix(sContourRegions) > 0)', .5, 'LineWidth', 1);

if (islogical(vnLabelRegions))
   vnLabelRegions = find(vnLabelRegions(:));
end

for (nRegion = reshape(vnLabelRegions, 1, []))
   [y, x] = ind2sub(vnStackSize, sCellRegions.PixelIdxList{nRegion});
   text(mean(y), mean(x), num2str(nRegion), ...
      'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'color', 'w');
end

PlotScaleBar(fPixPerUM, 25, 'tr', 'w-', 'LineWidth', 8);

if (nargout == 0)
   clear hFigure;
end


% --- END of MakeROIOverviewFigure.m ---
