function [hFigure,scont] = PlotROIs(sCellRegions, fPixPerUM, vnStackSize, vnLabelRegions)
% PlotROIs - FUNCTION
%
% Usage: [hFigure] = PlotRois(sCellRegions, fPixPerUM)
%
% sCellRegions is output from ReadImageJROI
%
% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 23rd May, 2011
% modified 120621 kampa@hifo.uzh.ch
%
% -- Defaults

if (~exist('vnLabelRegions', 'var'))
   vnLabelRegions = 1:sCellRegions.NumObjects;
end

if (~exist('vnStackSize', 'var'))
    vnStackSize = [256 64];
end

newplot;
hFigure = gcf;

sContourRegions = sCellRegions;
sContourRegions.NumObjects = nnz(vnLabelRegions);
sContourRegions.PixelIdxList = sContourRegions.PixelIdxList(vnLabelRegions);

for i=1:length(vnLabelRegions)
    [C,scont(i)] = contour((labelmatrix(sContourRegions) == i)', .5, 'LineWidth', 1,'Color','k');
    hold on;
end;
 hold off;
 
if (islogical(vnLabelRegions))
   vnLabelRegions = find(vnLabelRegions(:));
end

for (nRegion = reshape(vnLabelRegions, 1, []))
   [y, x] = ind2sub(vnStackSize, sCellRegions.PixelIdxList{nRegion});
   text(mean(y), mean(x), num2str(nRegion), ...
      'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'color', 'k');
end

PlotScaleBar(fPixPerUM, 5, 'tr', 'k-', 'LineWidth', 8);

set(gca,'YDir','reverse')

end