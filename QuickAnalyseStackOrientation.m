function [mbRespMask, fsStack] = QuickAnalyseStackOrientation(strBaseName, vnStackNumbers, mbUseMask, vfOrientations, fsStack)

% QuickAnalyseStackOrientation - FUNCTION 
%
% Usage: [mbRespMask, fsStack] = QuickAnalyseStackOrientation(strBaseName, vnStackNumbers, mbUseMask, vfOrientations, fsStack)

tic;

% -- Defaults

DEF_vfOrientations = linspace(0, pi, 9);
DEF_vfOrientations = DEF_vfOrientations(1:8);

nMinRegionSize = 9;

%% -- Check arguments

if (nargin < 2)
   disp('*** QuickAnalyseStackOrientation: Incorrect usage');
   help QuickAnalyseStackOrientation;
   return;
end

if (~exist('mbUseMask', 'var') || isempty(mbUseMask))
   mbUseMask = [];
end

if (~exist('vfOrientations', 'var') || isempty(vfOrientations))
   vfOrientations = DEF_vfOrientations;
end


%% -- Generate stack filenames

nNumFiles = numel(vnStackNumbers);
cstrStackFilenames = cell(nNumFiles, 1);

for (nFile = 1:nNumFiles) %#ok<FORPF>
   cstrStackFilenames{nFile} = sprintf('%s%03d.fcs', strBaseName, vnStackNumbers(nFile));
end

if (~exist('fsStack', 'var') || isempty(fsStack))
   % - Make a stack and align it
   fsStack = FocusStack(cstrStackFilenames);
   fprintf(1, '--- QuickAnalyseStackOrientation: Aligning stack...\n');
%    fsStack.Align(2, false, 10);
end

vnStackSize = size(fsStack);
% mbAlignMask = fsStack.GetAlignedMask;


%% -- Set up frame indices

tBlankTime = 10;        % sec
tUseBlank = 14;        % sec
tFrameRate = 3.91;      % Hz
tEndTrialBlankTime = 0; % sec

tInterStim = 5;      % sec
tStimTime = 5;       % sec
tUseStart = 0;       % sec
tUseEnd = 5;         % sec
nStimSkipStart = round(tUseStart * tFrameRate) + 1;
nStimUseEnd = round(tUseEnd * tFrameRate);

nNumStim = numel(vfOrientations);

nNumRepeats = nNumFiles;
nFramesPerTrial = vnStackSize(3) / nNumRepeats;
tTimeTrace = (0:nFramesPerTrial-1) ./ tFrameRate;

vnBlankFrames = round([1:(tUseBlank*tFrameRate) nFramesPerTrial-(tEndTrialBlankTime*tFrameRate):nFramesPerTrial]);
vnBlankFrames = repmat(vnBlankFrames, nNumRepeats, 1);
vnBlankFrames = vnBlankFrames + repmat(((0:nNumRepeats-1) * nFramesPerTrial)', 1, size(vnBlankFrames, 2));
vnBlankFrames = vnBlankFrames';
vnBlankFrames = vnBlankFrames(:);

for (nStim = nNumStim:-1:1)
   tStimStart = tBlankTime + (nStim)*tInterStim + (nStim-1)*tStimTime;
   tStimEnd = tStimStart + tStimTime;
   
   cvnStimFrames{nStim} = round((tStimStart*tFrameRate:tStimEnd*tFrameRate)); %#ok<*AGROW,*SAGROW>
   cvnStimFrames{nStim} = cvnStimFrames{nStim}(nStimSkipStart:nStimUseEnd);
   
   cvnStimFrames{nStim} = repmat(cvnStimFrames{nStim}, nNumRepeats, 1);
   cvnStimFrames{nStim} = cvnStimFrames{nStim} + repmat(((0:nNumRepeats-1) * nFramesPerTrial)', 1, size(cvnStimFrames{nStim}, 2));
   cvnStimFrames{nStim} = cvnStimFrames{nStim}';
   cvnStimFrames{nStim} = cvnStimFrames{nStim}(:);
end

nFramesPerStim = tStimTime * tFrameRate;


%% -- Find responsive pixels

fprintf(1, '--- QuickAnalyseStackOrientation: Locating responsive regions...\n');

sWarn = warning('off', 'FocusStack:UnalignedStack');

if (isempty(mbUseMask))
   mbUseMask = false(vnStackSize(1:2));
end

[mbRespMask] = ResponsiveMask(fsStack, vnBlankFrames, cvnStimFrames, 1, 1e-3);

% mbRespMask = mbAlignMask & (mbRespMask | mbUseMask);
mbRespMask = (mbRespMask | mbUseMask);

sRespRegions = bwconncomp(mbRespMask, 8);


%% -- Filter out small regions

vnRegionSizes = cellfun(@numel, sRespRegions.PixelIdxList);
vbAcceptRegion = vnRegionSizes >= nMinRegionSize;
sRespRegions.PixelIdxList = sRespRegions.PixelIdxList(vbAcceptRegion);
sRespRegions.NumObjects = numel(sRespRegions.PixelIdxList);


%% - Extract aligned traces for all identified responsive pixels

fprintf(1, '--- QuickAnalyseStackOrientation: Extracting region responses...\n');

[vfBlankMeans, nul, mfStimMeans, nul, mfRegionTraces] = ...
   ExtractRegionResponses(fsStack, sRespRegions, vnBlankFrames, cvnStimFrames, 1);

% - Average trials (other stims)
mfRegionTracesAvg = zeros(sRespRegions.NumObjects, nFramesPerTrial);

for (nTrial = 1:nNumRepeats)
   mfRegionTracesAvg(:, 1:nFramesPerTrial) = mfRegionTracesAvg(:, 1:nFramesPerTrial) + ...
      mfRegionTraces(:, (nTrial-1)*nFramesPerTrial+1:nTrial*nFramesPerTrial);
end

mfRegionTracesAvg = mfRegionTracesAvg(:, 1:nFramesPerTrial);

% - Compute average over trials, convert to dF/F
mfRegionTraces = mfRegionTraces ./ repmat(vfBlankMeans, 1, size(mfRegionTraces, 2));
mfRegionTracesAvg = mfRegionTracesAvg ./ nNumRepeats ./ repmat(vfBlankMeans, 1, size(mfRegionTracesAvg, 2));
mfStimMeans = mfStimMeans ./ repmat(vfBlankMeans, 1, size(mfStimMeans, 2));

% - Make a stack image
vnStackSize = size(fsStack);
mfGreenMean = fsStack.SummedAlignedFrames(:, :, :, 1) ./ vnStackSize(3);
mfRedMean = fsStack.SummedAlignedFrames(:, :, :, 2) ./ vnStackSize(3);
imResponse = cat(3, mfRedMean', mfGreenMean', zeros(vnStackSize(1:2)));


%% Make a nice plot showing signals

vhFigures = figure;

image(uint8(imResponse));
hold on;
contour((labelmatrix(sRespRegions) > 0)', .5, 'k-', 'LineWidth', 2);
title(sprintf('%s [%d:%d]', strBaseName, min(vnStackNumbers), max(vnStackNumbers)), 'FontSize', 24);

axis equal tight off;
set(gcf, 'Color', 'w');


for (nRegion = 1:sRespRegions.NumObjects)
   [y, x] = ind2sub(vnStackSize(1:2), sRespRegions.PixelIdxList{nRegion});
   text(mean(y), mean(x), num2str(nRegion), ...
      'FontSize', 12, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'w');
end

for (nRegion = 1:sRespRegions.NumObjects)
   if (mod(nRegion-1, 5) == 0)
      vhFigures(end+1) = figure; %#ok<AGROW>
      set(gcf, 'Color', 'w');
   end
   
   % - Create a subplot
   subplot(5, 1, mod(nRegion-1, 5)+1);
   
   fYMin = min(mfRegionTraces(nRegion, :));
   fYMax = max(mfRegionTraces(nRegion, :));
   vfYLims = [fYMin fYMax];
   
   % - Plot the stimulus times
   for (nStim = 1:nNumStim)
      tStimStart = tBlankTime + (nStim)*tInterStim + (nStim-1)*tStimTime;
      tStimEnd = tStimStart + tStimTime;
      vtStimTimes = [tStimStart tStimEnd];
      hPatch = fill(vtStimTimes([1 1 end end 1]), vfYLims([1 2 2 1 1]), [0.85 0.85 0.85]);
      set(hPatch, 'LineStyle', 'none');
      hold on;
   end
   
   for (nTrial = 1:nNumRepeats)
      plot(tTimeTrace, mfRegionTraces(nRegion, (nTrial-1)*nFramesPerTrial+1:nTrial*nFramesPerTrial), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
   end
   
   % - Plot the average fluorescence trace for the region
   plot(tTimeTrace, mfRegionTracesAvg(nRegion, :), 'r-', 'LineWidth', 2);
   
   % - Plot a dF/F scale bar
   fdFFScale = 0.1;
   vfBarLims = vfYLims(2) - [fdFFScale 0] - 0.05;
   plot(tTimeTrace([5 5]), vfBarLims, 'k-', 'LineWidth', 5);
   
   % - Set the plot properties
   axis([0 tTimeTrace(end) fYMin fYMax]);
   box on;
   
   if (nRegion ~= sRespRegions.NumObjects) && (mod(nRegion-1, 5) ~= 4)
      set(gca, 'XTick', []);
   end
   
   set(gca, 'Color', 'none', 'YTick', []);
   ylabel(gca, nRegion, 'FontSize', 14);
end


%% Assign a preferred orientation to each region

vfOriVectors = cos(vfOrientations*2) + 1i*sin(vfOrientations*2);

mfOriResponses = mfStimMeans .* repmat(vfOriVectors, sRespRegions.NumObjects, 1);
vfPrefOri = angle(mean(mfOriResponses, 2)) / 2;
vfOriSelStrength = abs(mean(mfOriResponses, 2));

vfPrefOri(vfPrefOri < 0) = vfPrefOri(vfPrefOri < 0) + pi;
vfPrefOri(vfPrefOri > pi) = vfPrefOri(vfPrefOri > pi) - pi;


%% Make a plot showing distribution of responses

figure;
hist(vfPrefOri * 180/pi);


%% -- Show some timing feedback

warning(sWarn);

fprintf(1, '--- QuickAnalyseStackOrientation: Analysed [%d] frames in [%d] seconds.\n', vnStackSize(3), round(toc));

% --- END of QuickAnalyseStackOrientation.m ---
