function [fsData, sQuickAnalysis, vhFigures] = QuickAnalyseStack( strBaseFilename, vnFileNumbers, ...
                                                                  bAlign, bIdentifyBlack, bSaveFigures)
                                           
% QuickAnalyseStack - FUNCTION Perform a quick analysis of a focus stack
%
% Usage: [fsData, sQuickAnalysis, vhFigures] = QuickAnalyseStack( strBaseFilename, vnFileNumbers, ...
%                                                                 bAlign, bIdentifyBlack, bSaveFigures)

if (nargin < 2)
   disp('*** QuickAnalyseStack: Incorrect usage.');
   help QuickAnalyseStack;
   return;
end

if (~exist('bAlign', 'var'))
   bAlign = false;
end

if (~exist('bIdentifyBlack', 'var'))
   bIdentifyBlack = false;
end

if (~exist('bSaveFigures', 'var'))
   bSaveFigures = false;
end

if (isa(strBaseFilename, 'FocusStack'))
   fsData = strBaseFilename;
   vnStackSize = size(fsData);
   
else
   %% -- Build the stack
   disp('--- QuickAnalyseStack: Building stack...');
   cstrFilenames = {};
   for (nFile = vnFileNumbers)
      cstrFilenames{end+1} = sprintf('%s%03d.fcs', strBaseFilename, nFile); %#ok<AGROW,SAGROW>
   end
   
   % - Make a focus stack
   fsData = FocusStack(cstrFilenames);
   vnStackSize = size(fsData);
end

if (bAlign && isempty(fsData.mfFrameShifts))
   fsData.Align(2, false, 10);
end

if (bIdentifyBlack)
   disp('--- QuickAnalyseStack: Assign a black region...');
   fsData.DefineBlackRegion();
end

% - Turn off warnings
w = warning('off', 'FocusStack:IncompleteInformation');

if (~bAlign)
   warning('off', 'FocusStack:UnalignedStack');
end


%% -- Locate neurons

disp('--- QuickAnalyseStack: Locating cells (ROIs)...');

[sCellRegions] = FindCells_GRChannels(fsData);


%% -- Extract average response

disp('--- QuickAnalyseStack: Extracting average response...');
fsData.BlankNormalisation('none');
tfAvgSignal = fsData.SummedAlignedFrames(:, :, :, :) ./ vnStackSize(3);


%% -- Plot a figure of the responsive regions

vnStackSize = size(fsData);

imAvg = tfAvgSignal;

% - Transpose image and normalise
imAvg(:, :, 1) = tfAvgSignal(:, :, 2)' - min(min(tfAvgSignal(:, :, 2)));
imAvg(:, :, 1) = imAvg(:, :, 1) ./ max(max(imAvg(:, :, 1)));
imAvg(:, :, 2) = tfAvgSignal(:, :, 1)' - min(min(tfAvgSignal(:, :, 1)));
imAvg(:, :, 2) = imAvg(:, :, 2) ./ max(max(imAvg(:, :, 2)));
imAvg(:, :, 3) = 0;

hFigure1 = figure;
subplot(1, 2, 1);
image(imAvg);
axis equal tight off;

subplot(1, 2, 2);
image(imAvg);
axis equal tight off;
hold on;
contour((labelmatrix(sCellRegions) > 0)', .5, 'LineWidth', 1);
% alpha(hLabels, double(mbResponsive')*0.2);

for (nRegion = 1:sCellRegions.NumObjects)
   [y, x] = ind2sub(vnStackSize(1:2), sCellRegions.PixelIdxList{nRegion});
   text(mean(y), mean(x), num2str(nRegion), ...
      'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'color', 'w');
end
drawnow;


%% -- Extract region responses

disp('--- QuickAnalyseStack: Extracting region responses...');

[nul, nul, nul, nul, mfRegionTraces] = ...
   ExtractRegionResponses(fsData, sCellRegions, 1);


%% -- Extract neuropil response

disp('--- QuickAnalyseStack: Extracting neuropil response...');

mbNeuropilMask = ~labelmatrix(sCellRegions);
vfNeuropilTrace = fsData.MeanPixels(mbNeuropilMask, :, 1);

sTemp = bwconncomp(mbNeuropilMask);
sNeuropilRegion = sTemp;
sNeuropilRegion.NumObjects = 1;
sNeuropilRegion.PixelIdxList = {vertcat(sTemp.PixelIdxList{:})};
hFigure2 = PlotRegionResponses(vfNeuropilTrace, fsData, sNeuropilRegion, [], 0.1, 5);
title('Neuropil mean response');


%% -- Plot responses

vhFigures = PlotRegionResponses(mfRegionTraces, fsData, sCellRegions, [], 0.1, 5, 30);

% - Restore warnings
warning(w);

sQuickAnalysis.sCellRegions = sCellRegions;
sQuickAnalysis.mfRegionTraces = mfRegionTraces;
sQuickAnalysis.vfNeuropilTrace = vfNeuropilTrace;

vhFigures = [hFigure1 hFigure2 vhFigures];


% return;

%% Write figures to PDF

if (bSaveFigures)
   strOldDir = cd('../analysis');
   
   strAnalysisName = [strBaseFilename 'analysis_' num2str(vnFileNumbers(1)) '-' num2str(vnFileNumbers(end))];
   
   for (hFigure = vhFigures)
      export_to_pdf(hFigure, [strAnalysisName '_' num2str(hFigure) '.pdf']);
   end
   
   save([strAnalysisName '_quick.mat'], 'fsData', 'sQuickAnalysis');
   
   cd(strOldDir);
end


% --- END of QuickAnalyseStack.m ---

