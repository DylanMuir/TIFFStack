function [fsData, sAnalysis, vhFigures] = AnalyseStack( strBaseFilename, vnFileNumbers, ...
                                             bAlign, vtStimDurations, fAlpha, ...
                                             mtUseStimDurations, mtBlankDurations, mtStimLabelTimes, tTimeShift, ...
                                             cvnSequenceIDs, fPixelsPerUM, bPlotFigures)
                                          
% AnalyseStack - FUNCTION Perform a series of analyses on a Focus stack
%
% Usage: [fsData, sAnalysis, vhFigures] = AnalyseStack( strBaseFilename, vnFileNumbers, ...
%                                            bAlign, vtStimDurations, <fAlpha>, ...
%                                            <mtUseStimDurations, mtBlankDurations, mtStimLabelTimes, tTimeShift>, ...
%                                            <cvnSequenceIDs, fPixelsPerUM, bPlotFigures>)
%

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 13th January, 2011

if (~exist('bPlotFigures', 'var') || isempty(bPlotFigures))
   bPlotFigures = true;
end

if (isa(strBaseFilename, 'FocusStack'))
   fsData = strBaseFilename;
   vnStackSize = size(fsData);
   
else
   %% -- Build the stack
   cstrFilenames = {};
   for (nFile = vnFileNumbers)
      cstrFilenames{end+1} = sprintf('%s%03d.fcs', strBaseFilename, nFile); %#ok<AGROW,SAGROW>
   end
   
   % - Make a focus stack
   fsData = FocusStack(cstrFilenames);
   vnStackSize = size(fsData);
   
  
   %% -- Embed stimulus sequence, if necessary
   
   if (exist('cvnSequenceIDs', 'var') && ~isempty('cvnSequenceIDs'))
      cvnStackSequenceIDs = fsData.cvnSequenceIDs;
      vbAssignSequence = cellfun(@(c)(isnan(c(1))), cvnStackSequenceIDs);
      
      cvnStackSequenceIDs(vbAssignSequence) = cvnSequenceIDs(vbAssignSequence);
      
      if (any(vbAssignSequence))
         % - Assign stimulus sequence
         fsData.cvnSequenceIDs = cvnStackSequenceIDs;
      end
   end
   
   % - Assign stimulus durations
   fsData.vtStimulusDurations = vtStimDurations;

   if (~exist('mtUseStimDurations', 'var') || isempty(mtUseStimDurations))
      mtUseStimDurations = fsData.mtStimulusUseTimes;%[zeros(numel(vtStimDurations, 1)) reshape(vtStimDurations, [], 1)];
   else
      fsData.mtStimulusUseTimes = mtUseStimDurations;
   end
   
   if (~exist('mtStimLabelTimes', 'var'))
      mtStimLabelTimes = [];
   end
   
   if (exist('tTimeShift', 'var') && ~isempty(tTimeShift))
      % - Shift stim start times
      vtStimStarts = fsData.vtStimulusStartTimes;
      fsData.vtStimulusStartTimes = vtStimStarts + tTimeShift;
      
      % - Shift stim use times and plot times
      mtUseStimDurations = mtUseStimDurations - tTimeShift;
      mtUseStimDurations(mtUseStimDurations < 0) = 0;
      mtBlankDurations = mtBlankDurations - tTimeShift;
      mtBlankDurations(mtBlankDurations < 0) = 0;
      mtStimLabelTimes = mtStimLabelTimes - tTimeShift;
      mtStimLabelTimes(mtStimLabelTimes < 0) = 0;
      
      fsData.mtStimulusUseTimes = mtUseStimDurations;
   end
   
   if (exist('fPixelsPerUM', 'var') && ~isempty(fPixelsPerUM))
      fsData.fPixelsPerUM = fPixelsPerUM;
   end
   
   
   %% -- Align the stack
   
   if (bAlign && isempty(fsData.mfFrameShifts))
      fsData.Align(2, false, 10);
   end
   
   
   %% -- Define a black region
   
   disp('--- AnalyseStack: Assign a black region...');
   fsData.DefineBlackRegion();   
end


%% -- Extract average response

disp('--- AnalyseStack: Extracting average response...');
fsData.BlankNormalisation('none');
tfAvgSignal = fsData.SummedAlignedFrames(:, :, :, :) ./ vnStackSize(3);


%% -- Get stack frame info and assign blank frames

disp('--- AnalyseStack: Extracting blank frames...');

% - Get an alignment mask
mbAlignMask = fsData.GetAlignedMask;

[vtGlobalTime, ...
   vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
   vnStimulusSeqID, vtTimeInStimPresentation] = ...
   FrameStimulusInfo(fsData, 1:size(fsData, 3));

% - Turn off blank warning and unaligned warning
wOld = warning('off', 'FocusStack:FilteringBlankFrames');
warning('off', 'FocusStack:UnalignedStack');

% - Turn off DFF conversion, just in case
fsData.BlankNormalisation('none');

for (nBlock = 1:numel(fsData.cstrFilenames))
   % - Get the average blank frame for this block
   vbBlockFrames = vnBlockIndex == nBlock;
   vbBlockBlankStimFrames = vbBlockFrames & (vnStimulusSeqID == 1);
   vbBlockBlankFrames = vbBlockBlankStimFrames & ...
      (vtTimeInStimPresentation >= mtUseStimDurations(1, 1)) & ...
      (vtTimeInStimPresentation <= mtUseStimDurations(1, 2));
   
   tfBlankFrames = double(fsData.AlignedStack(:, :, vbBlockBlankFrames, 1));
   mfThisBlankAvg = mean(tfBlankFrames, 3);
   mfThisBlankStd = std(tfBlankFrames, [], 3);
   
   % - Assign the blank frame to the whole block (by default)
   fsData.AssignBlankFrame(cat(3, mfThisBlankAvg, mfThisBlankStd), vbBlockFrames, 1);
   
   % - Assign separate blank frames for each other stimulus segment
   for (nStimSeqID = 2:numel(vtStimDurations))
      % - Find all frames for this presentation
      vbPresFrames = vbBlockFrames & (vnStimulusSeqID == nStimSeqID);
      
      % - Find blank frames for this presentation
      vbPresBlankFrames =  vbPresFrames & ...
                           (vtTimeInStimPresentation >= mtBlankDurations(nStimSeqID, 1)) & ...
                           (vtTimeInStimPresentation <= mtBlankDurations(nStimSeqID, 2));

      % - Extract these blank frames
      tfPresBlank = fsData.AlignedStack(:, :, vbPresBlankFrames, 1);
      
      % - Assign only the mean; Std.Dev is taken from block blank
      fsData.AssignBlankFrame(cat(3, mean(tfPresBlank, 3), mfThisBlankStd), vbPresFrames, 1);
   end
end

% - Restore warnings
warning(wOld);



%% -- Classify responsive regions

disp('--- AnalyseStack: Locating cells...');
[sCellRegions] = FindCells_GRChannels(fsData);

if (~exist('fAlpha', 'var'))
   fAlpha = 0.05;
end

if (fsData.nNumStimuli == 1)
   vnUseStimIDs = 1;
else
   vnUseStimIDs = 2:fsData.nNumStimuli;
end


disp('--- AnalyseStack: Classifying responsive regions...');
[sRespRegions, nul, mfMaxZScores, mfMaxFiltZScores, tfStimZFilt] = ...
   ResponsiveMask(fsData, vnUseStimIDs, 1, fAlpha, 6, 200);

% - Find max Z score for each region
vfRegionZScore = cellfun(@(r)(max(mfMaxFiltZScores(r))), sCellRegions.PixelIdxList);

% - Convert to alphas
vfRegionAlpha = 1-normcdf(vfRegionZScore);

% - Estimate a neuropil-rejecting Z score threshold
[fResponsiveMaxZThresh, fResponsiveMaxEquivAlpha] = ...
   EstimateResponsivenessThreshold(mfMaxFiltZScores, sCellRegions, fAlpha);

[fResponsiveAnyZThresh, fResponsiveAnyEquivAlpha] = ...
   EstimateResponsivenessThreshold(tfStimZFilt, sCellRegions, fAlpha);


%% -- Extract aligned traces for all identified responsive pixels

disp('--- AnalyseStack: Extracting region responses...');

fsData.BlankNormalisation('subtractive');

[vfBlankMeans, vfBlankStds, mfStimMeans, mfStimStds, mfRegionTraces, mnStimSampleSizes] = ...
   ExtractRegionResponses(fsData, sCellRegions, 1);


%% -- Re-order by max stimulus Z score

mfStimZScores = mfStimMeans ./ (repmat(vfBlankStds', 1, size(mfStimMeans, 2)) ./ sqrt(mnStimSampleSizes));

[vfSortedZScores, vnSort] = sort(max(mfStimZScores, [], 2), 'descend');

sCellRegions.PixelIdxList = sCellRegions.PixelIdxList(vnSort);
sCellRegions.vfRegionZScore = vfSortedZScores;
sCellRegions.vfRegionAlpha = 1 - normcdf(vfSortedZScores);

vfBlankMeans = vfBlankMeans(vnSort);
vfBlankStds = vfBlankStds(vnSort);
mfStimMeans = mfStimMeans(vnSort, :);
mfStimStds = mfStimStds(vnSort, :);
mfRegionTraces = mfRegionTraces(vnSort, :);
mnStimSampleSizes = mnStimSampleSizes(vnSort, :);
mfStimZScores = mfStimZScores(vnSort, :);

%% -- Plot figures

% - Make a mask matrix
mbResponseMask = (labelmatrix(sCellRegions) > 0);
   
if (bPlotFigures)

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
   image(imAvg);
   axis equal tight off;
   
   hFigure2 = figure;
   image(imAvg);
   axis equal tight off;
   hold on;
   contour((labelmatrix(sCellRegions) > 0)', .5, 'LineWidth', 1);
   
   for (nRegion = 1:sCellRegions.NumObjects)
      [y, x] = ind2sub(vnStackSize(1:2), sCellRegions.PixelIdxList{nRegion});
      text(mean(y), mean(x), num2str(nRegion), ...
         'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'color', 'w');
   end
   
   
   %% -- Make a plot of the response traces for different stimuli
   
   vhFigures = PlotRegionResponses(mfRegionTraces, fsData, sCellRegions, mtStimLabelTimes, 0.1, 5, 100);
   vhFigures = [hFigure1 hFigure2 vhFigures];
else
   vhFigures = [];
end

%% -- Store analysis data

sAnalysis.sCellRegions = sCellRegions;
sAnalysis.sRespRegions = sRespRegions;
sAnalysis.mbResponseMask = mbResponseMask;
sAnalysis.mfMaxZScores = mfMaxZScores;
sAnalysis.mfMaxFiltZScores = mfMaxFiltZScores;
sAnalysis.tfStimZFilt = tfStimZFilt;
sAnalysis.fResponsiveMaxZThresh = fResponsiveMaxZThresh;
sAnalysis.fResponsiveMaxEquivAlpha = fResponsiveMaxEquivAlpha;
sAnalysis.vfBlankMeans = vfBlankMeans;
sAnalysis.vfBlankStds = vfBlankStds;
sAnalysis.mfStimMeans = mfStimMeans;
sAnalysis.mfStimStds = mfStimStds;
sAnalysis.mfStimZScores = mfStimZScores;
sAnalysis.mnStimSampleSizes = mnStimSampleSizes;


%% -- Save analysis data

if (ischar(strBaseFilename))
   save(sprintf('%s%d-%d_analysis.mat', strBaseFilename, min(vnFileNumbers), max(vnFileNumbers)), ...
      'fsData', 'sAnalysis');
end

% --- END of AnalyseStack.m ---
