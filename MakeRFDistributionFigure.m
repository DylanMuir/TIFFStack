function hFigure = MakeRFDistributionFigure(sRFAnalysis, vfStimOffsetDeg, fStimRadiusDeg, vnPlotRegions)

% MakeRFDistributionFigure - FUNCTION Make a figure showing the location and extent of estimated RFs
%
% Usage: hFigure = MakeRFDistributionFigure(sRFAnalysis, vfStimOffsetDeg, fStimRadiusDeg, vnPlotRegions)
%
% 

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 25th March, 2013

if (~exist('vnPlotRegions', 'var') || isempty(vnPlotRegions))
   vnPlotRegions = 1:numel(sRFAnalysis.vbReliableROIs);
end

if (islogical(vnPlotRegions))
   vnPlotRegions = find(vnPlotRegions);
end

hFigure = figure;

% - Plot an outline of the display screen
% plot(vfStimOffsetDeg(1), vfStimOffsetDeg(2), 'kx', 'MarkerSize', 8, 'LineWidth', 2);
% hold on;
% plot(sRFAnalysis.vfScreenSizeDeg(1) * [-0.5 0.5 0.5 -0.5 -0.5], sRFAnalysis.vfScreenSizeDeg(2) * [-0.5 -0.5 0.5 0.5 -0.5], 'g:', 'LineWidth', 2);

% - Make an RF estimate sample mesh
vfX = linspace(-sRFAnalysis.vfScreenSizeDeg(1)/2, sRFAnalysis.vfScreenSizeDeg(1)/2, ...
   size(sRFAnalysis.vsRFEstimates(1).mfRFImage, 1));
vfY = linspace(-sRFAnalysis.vfScreenSizeDeg(2)/2, sRFAnalysis.vfScreenSizeDeg(2)/2, ...
   size(sRFAnalysis.vsRFEstimates(1).mfRFImage, 2));
[mfX, mfY] = ndgrid(vfX, vfY);

w=warning('off', 'MATLAB:contour:ConstantData');

% - Which ROIs should be accepted
vbAcceptROIs = sRFAnalysis.vbReliableROIs & sRFAnalysis.vbResponsiveROIs & sRFAnalysis.vbAcceptRFLocation;

% - Plot contours for each RF
for (nRegion = vnPlotRegions(:)')
   if (vbAcceptROIs(nRegion))
      strLinespec = 'k-';
   else
      strLinespec = 'r--';
   end
   
   fMaxLevel = max(sRFAnalysis.vsRFEstimates(nRegion).mfRFImage(:));
   
   contour(mfX, mfY, sRFAnalysis.vsRFEstimates(nRegion).mfRFImage, fMaxLevel * [.5 .5], strLinespec, 'LineWidth', 1);
   hold on;
end
      
   
% - Draw a circle for reference
if (exist('fStimRadiusDeg', 'var') && ~isempty(fStimRadiusDeg))
   vfAngle = linspace(0, 2*pi, 200);
   vfX = vfStimOffsetDeg(1) + fStimRadiusDeg * cos(vfAngle);
   vfY = vfStimOffsetDeg(2) + fStimRadiusDeg * sin(vfAngle);
   plot(vfX, vfY, 'k:', 'LineWidth', 2);
end

axis equal tight

% - Restore warnings
warning(w);

% --- END of MakeRFDistributionFigure.m ---
