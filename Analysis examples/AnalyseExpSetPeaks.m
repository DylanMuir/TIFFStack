function [sAnalysis, vhFigures] = AnalyseExpSetPeaks( sExpSpec, bAlign, bAssignBlack, vnBlackPixels, bPlotFigures, ...
                                                      vbForceAnalysisSections, sManualRegions, mfAlignmentImage, ...
                                                      tBaseTimeShift, fForcePixelsPerUM, cvnForceStackSequenceIDs, tPeakSearchTime)

% AnalyseExpSetPeaks - FUNCTION Extract calcium responses from raw 2P imaging data, using the peak response
%
% Usage: [sAnalysis, vhFigures] = AnalyseExpSetPeaks( sExpSpec, < bAlign, bAssignBlack, vnBlackPixels, bPlotFigures, ...
%                                                     vbForceAnalysisSections, sManualRegions, mfAlignmentImage, ...
%                                                     tBaseTimeShift, fForcePixelsPerUM, cvnForceStackSequenceIDs, tPeakSearchTime>)
%
% This function extracts calcium responses from a set of raw two-photon
% imaging stacks, in response to a sequence of stimuli. The response
% measure is dF/F, measured as the peak calcium deflection following
% stimulus onset. Baseline F is defined as the median response during a
% pre-stimulus baseline period.
%
% 'sExpSpec' is a structure that defines where the 2P raw data is located,
% and other experimental parameters. It has the format:
%
% sExpSpec.strExpName - String containing the name of this experiment. Free text.
%         .strExpDate - String containing the date of this experiment. Free text.
%         .strBasePath - File path to the base directory containing experimental data and analysis results
%         .strExpDataLoc - Subdirectory under base path that contains data for this experiment
%         .strDataSubdir - Subdirectory under data loc that contains raw data
%         .strAnalysisSubdir - Subdirectory under base path to store analysis results
%         .vnUseFileIndices - Vector of data file indices to use for analysis
%         .strROISetFilename - ImageJ ROISet File under data subdir that defines ROIs to analyse
%         .strFileType - String containg the file extension of the raw data ('.fcs', '.bin', '.tif', etc)
%         .nNumTrueStim - Integer number of non-blank stimuli
%         .nBlankStimID - Index of blank stimulus
%         .tStimDuration - Duration of non-blank stimuli, in seconds
%         .tBlankStimDuration - Duration of blank stimulus, in seconds
%
% The raw data files are assumed to be named 'strExpDate'_'index''strFileType', 
% where 'index' is one of the files defined by 'vnUseFileIndices'.
%
% 'bAlign' is a boolean flag indicating whether or not the stack should be
% aligned before analysis. By default, AlignTrials is used with either
% 'mfAlignmentImage' as a reference, or with a reference window of frames
% 10:50, and with 10-times subpixel alignment.
%
% 'bAssignBlack' is a boolean flag indicating whether or not the black
% trace should be extracted from the stack before analysis.
% 'vnBlackPixels' is a vector list of linear pixel indices that define
% which pixels should be used to extract the black trace for the stack.
%
% 'bPlotFigures' is a boolean flag determining whether or not analysis
% figures should be produced and saved to disk.
%
% 'vbForceAnalysisSections' is a vector of boolean flags that determine
% which elements of the analysis should not be taken from cached results,
% but should be recalculated:
%  [bStartFromScratch
%   bLoadOrIdentifyROIs
%   bExtractRegionResponses
%   bExtractAverageResponseImage
%   bEstimateResponseThresholds]
%   
% 'sManualRegions' is a matlab regions structure that overrides the ROIs
% defined in 'sExpSpec'
%
% 'mfAlignmentImage' is an image used as a reference for aligning the
% stack. This argument can also be supplied as 'fhAlignmentFunction': a
% function handle that is called to perform a custom alignment step. For
% example,
%    mfAlignmentImage = @(fs)fs.Align(1, false, 10, mfRefImage, 1, [1/200 1/20], 7)
%
% 'tBaseTimeShift' - Time in seconds to shift analysis time relative to
% true stack time. Positive values mean that frames used for analysis are
% taken from later in imaging time.
%
% 'fForcePixelsPerUM' is a value for the stack calibration, if that
% information is incorrect or not provided by the stack. The value is in
% pixels per micrometre.
%
% 'cvnForceStackSequenceIDs' is a cell array, each cell corresponding to
% one trial and one raw data source file. Each cell contains a vector that
% defines the order of stimuli presented during that trial. This argument
% is optional, and should be provided only if the data is not already
% present in the stack. Each cell should end in NAN, which means "skip to
% the next trial".
%
% 'tPeakSearchTime' is a time in seconds to search for the response peak,
% measured from stimulus onset.

% Author: Dylan Muir <muir@hifo.uzh.ch>, Morgane Roth <roth@hifo.uzh.ch>
% Created: 2013


% -- Check arguments

if (~exist('bAlign', 'var'))
   bAlign = true;
end

if (~exist('bAssignBlack', 'var'))
   bAssignBlack = true;
end

if (~exist('vnBlackPixels', 'var'))
   vnBlackPixels = [];
end

if (~exist('bPlotFigures', 'var'))
   bPlotFigures = true;
end

if (~exist('vbForceAnalysisSections', 'var'))
   vbForceAnalysisSections = false(1, 6);
end

if (~exist('tBaseTimeShift', 'var'))
   tBaseTimeShift = [];
end

if (~exist('cvnForceStackSequenceIDs', 'var'))
   cvnForceStackSequenceIDs = [];
end

if (~exist('tPeakSearchTime', 'var') || isempty(tPeakSearchTime))
   tPeakSearchTime = [];
end


%% -- Set up parameters

fAlpha = 0.05;
fBaselineProportion = 0.5;
fUsePreStimBlankProportion = 0.5;

% - Get analysis .mat filename
strAnalysisFile = GetExperimentAnalysisFilename(sExpSpec);

nNumOriStim = sExpSpec.nNumTrueStim + ~isempty(sExpSpec.nBlankStimID);

%% - Do we need to start from scratch?
if (~exist(strAnalysisFile, 'file') || vbForceAnalysisSections(1))
   % - Make a focus stack
   disp('--- AnalyseExpSetPeaks: Creating stack...');
   fsStack = GetExperimentFocusStack(sExpSpec);
   tPreStimBlank = fsStack.tBlankTime;
   
   if (isempty(tPreStimBlank))
      tPreStimBlank = 0;
   end
   
   vtStimulusDurations = repmat(sExpSpec.tStimDuration + tPreStimBlank, nNumOriStim, 1);
   vtStimulusDurations(sExpSpec.nBlankStimID) = sExpSpec.tBlankStimDuration;
   
   mtStimulusUseTimes = repmat(tPreStimBlank + [0 sExpSpec.tStimDuration], nNumOriStim, 1);
   mtStimLabelTimes = repmat(tPreStimBlank + [0 sExpSpec.tStimDuration], nNumOriStim, 1);
   mtStimLabelTimes(sExpSpec.nBlankStimID, :) = nan;
   mtBlankTimes = repmat(tPreStimBlank * [(1-fUsePreStimBlankProportion) 1], nNumOriStim, 1);
   if (~isempty(sExpSpec.nBlankStimID))
      mtBlankTimes(sExpSpec.nBlankStimID, :) = sExpSpec.tBlankStimDuration * [(1-fUsePreStimBlankProportion) 1];
      mtStimulusUseTimes(sExpSpec.nBlankStimID, :) = [0 sExpSpec.tBlankStimDuration];
   end
   vnDefSequenceIDs = [1:nNumOriStim nan];
   vnUseStimIDs = 1:nNumOriStim;
   
   if (~isempty(sExpSpec.nBlankStimID))
      vnUseStimIDs = vnUseStimIDs(vnUseStimIDs ~= sExpSpec.nBlankStimID);
   end
   
   
   % - Assign stimulus sequence, if required
if (exist('cvnForceStackSequenceIDs', 'var') && ~isempty(cvnForceStackSequenceIDs))
   if (~iscell(cvnForceStackSequenceIDs))
      cvnForceStackSequenceIDs = repmat({cvnForceStackSequenceIDs}, numel(fsStack.cstrFilenames), 1);
   end
   
   fsStack.cvnSequenceIDs = cvnForceStackSequenceIDs;
else
   % - Annotate stack
   cvnStackSequenceIDs = fsStack.cvnSequenceIDs;
   vbAssignSequence = cellfun(@(c)(isnan(c(1))), cvnStackSequenceIDs);
   
   cvnSequenceIDs = repmat({vnDefSequenceIDs}, 1, numel(fsStack.cstrFilenames));
   cvnStackSequenceIDs(vbAssignSequence) = cvnSequenceIDs(vbAssignSequence);
   
   if (exist('fForcePixelsPerUM', 'var'))
      fsStack.fPixelsPerUM = fForcePixelsPerUM;
   end
   
   if (any(vbAssignSequence))
      % - Assign stimulus sequence
      fsStack.cvnSequenceIDs = cvnStackSequenceIDs;
   end
end   
   fsStack.vtStimulusDurations = vtStimulusDurations;
   fsStack.mtStimulusUseTimes = mtStimulusUseTimes;
   
   % - Align stack
   if (bAlign)
      disp('--- AnalyseExpSetPeaks: Aligning stack...');
      if (exist('mfAlignmentImage', 'var'))
         if (isa(mfAlignmentImage, 'function_handle'))
            mfAlignmentImage(fsStack);
         else
            fsStack.AlignTrials(2, 10, mfAlignmentImage);
         end
      else
         fsStack.AlignTrials(2, 10, 10:50);
      end
   end
   
   % - Assign black trace
   if (bAssignBlack)
      fsStack.DefineBlackRegion(vnBlackPixels);
   end
   
   % - Turn off blank warning and unaligned warning
   wOld = warning('off', 'FocusStack:FilteringBlankFrames');
   warning('off', 'FocusStack:BlankContainsNaN');
   warning('off', 'FocusStack:UnalignedStack');
   
   % - Turn off DFF conversion, just in case
   fsStack.BlankNormalisation('none');
   
   % - Assign baseline
   fsStack = AssignBaselineDistribution(fsStack, vnUseStimIDs, sExpSpec.nBlankStimID, tBaseTimeShift, fBaselineProportion, mtBlankTimes);
   
   % - Restore warnings
   warning(wOld);
   
   % - Save stack
   save(strAnalysisFile, 'fsStack', 'mtStimLabelTimes', 'mtBlankTimes', 'vnUseStimIDs');
   
else
   % - Load analysis file
   disp('--- AnalyseExpSetPeaks: Loading analysis results...');
   load(strAnalysisFile);
end

%% Set up extraction function

% - Three points, peak and two following, (default) 2 sec search window
fhExtractionFun = ExtractPeak([0 1 2], 1, true, ceil(tPeakSearchTime ./ fsStack.tFrameDuration));


%% -- Get ROIs

if (~exist('sRegions', 'var') || vbForceAnalysisSections(2))
   disp('--- AnalyseExpSetPeaks: Loading/identifying ROIs...');
   
   if (exist('sManualRegions', 'var') && ~isempty(sManualRegions))
      % - Use manually-provided regions
      sRegions = sManualRegions;
      save(strAnalysisFile, 'sRegions', '-append');
      
      % - Force analysis to continue
      clear vfBlankStds;
      
   elseif (isfield(sExpSpec, 'strROISetFilename') && ~isempty(sExpSpec.strROISetFilename))
      % - Load ROIs
      cvsROIs = ReadImageJROI(fullfile(sExpSpec.strBasePath, sExpSpec.strExpDataLoc, sExpSpec.strROISetFilename));
      
      % - Make a regions structure
      sRegions = ROIs2Regions(cvsROIs, size(fsStack, 1:2));
      
      % - Save ROIs and regions
      save(strAnalysisFile, 'cvsROIs', 'sRegions', '-append');
      
   else
      % - Automatically identify ROIs
      disp('--- AnalyseExpSetPeaks: Identifying ROIs...');
      sRegions = FindCells_GRChannels(fsStack);
      
      save(strAnalysisFile, 'sRegions', '-append');
   end
   
   % - Reject ROIs that intersect with the edge or alignment mask
   mbAlignMask = fsStack.GetAlignedMask;
   mbEdgeRejectMask = ~mbAlignMask;
   mbEdgeRejectMask(1, :) = true;
   mbEdgeRejectMask(:, 1) = true;
   mbEdgeRejectMask(end, :) = true;
   mbEdgeRejectMask(:, end) = true;
   vnRejectPixels = find(mbEdgeRejectMask);
   vbRejectROI = cellfun(@(c)any(ismember(c, vnRejectPixels)), sRegions.PixelIdxList); %#ok<NASGU>

   save(strAnalysisFile, 'vbRejectROI', '-append');

   % - Force analysis to continue
   clear vfBlankStds;
end



%% -- Extract region responses

if (~exist('vfBlankStds', 'var') || vbForceAnalysisSections(3))
   disp('--- AnalyseExpSetPeaks: Extracting responses...');
   % - Read blank frames and data for cell regions
   fsStack.BlankNormalisation('none');
   [vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, tfTrialResponses, tnTrialSampleSizes, cfTrialTraces] = ...
      ExtractRegionResponses(fsStack, sRegions, sExpSpec.nBlankStimID, fhExtractionFun, tBaseTimeShift); %#ok<NASGU>
   
   % - Save results
   save(strAnalysisFile, 'vfBlankStds', 'mfStimMeanResponses', ...
      'mfStimStds', 'mfRegionTraces', 'tfTrialResponses', 'tnTrialSampleSizes', ...
      'cfTrialTraces', '-append');
   
   % - Clear Z-score ordering to force reordering
   clear mfStimZScores;
end


%% -- Extract average response

if (~exist('tfAvgSignal', 'var') || vbForceAnalysisSections(4))
   disp('--- AnalyseExpSetPeaks: Extracting average response...');
   fsStack.BlankNormalisation('none');
   tfAvgSignal = fsStack.SummedAlignedFrames(:, :, :, :) ./ size(fsStack, 3);
   
   % - Save results
   save(strAnalysisFile, 'tfAvgSignal', '-append');
end


%% -- Re-order by max stimulus Z score, detect reliable responses

if (~exist('mfStimZScores', 'var') || vbForceAnalysisSections(5))
   disp('--- AnalyseExpSetPeaks: Estimating response thresholds and Z-scores...');
   % - Estimate neuropil-rejecting Z score thresholds
   
   if (~isempty(sExpSpec.nBlankStimID))
      [sNeuropilZThresh] = ...
         ThresholdNeuropilResponse(fsStack, sRegions, vnUseStimIDs, fAlpha, sExpSpec.nBlankStimID, fhExtractionFun, tBaseTimeShift);
   else
      [vtGlobalTime, ...
         vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
         vnStimulusSeqID, vtTimeInStimPresentation, ...
         vnPresentationIndex, vbUseFrame] = ...
         FrameStimulusInfo(fsStack, 1:size(fsStack, 3), tBaseTimeShift);
      
      vbAllPresBlankFrames = false(size(vtGlobalTime));
      for (nStimSeqID = vnUseStimIDs)
         vbAllPresBlankFrames = vbAllPresBlankFrames | ((vnStimulusSeqID == nStimSeqID) & ...
            (vtTimeInStimPresentation >= mtBlankTimes(nStimSeqID, 1)) & ...
            (vtTimeInStimPresentation <= mtBlankTimes(nStimSeqID, 2)));
      end
      
      [sNeuropilZThresh] = ...
         ThresholdNeuropilResponse(fsStack, sRegions, vnUseStimIDs, fAlpha, vbAllPresBlankFrames, fhExtractionFun, tBaseTimeShift);
   end
      
   % - Compute Z-score against blank STDs
   mnStimSampleSizes = sum(tnTrialSampleSizes, 3);
   mfBlankStdsCorr = bsxfun(@rdivide, vfBlankStds', sqrt(mnStimSampleSizes));
   mfStimZScores = mfStimMeanResponses ./ mfBlankStdsCorr;
   
   tfBlankStdsCorr = bsxfun(@rdivide, vfBlankStds', sqrt(tnTrialSampleSizes));
   tfStimZScoresTrials = tfTrialResponses ./ tfBlankStdsCorr;
   
   % - Measure reliability
   tbSignificantTrial = tfStimZScoresTrials > sNeuropilZThresh.fTrialStimResp;
   %tbSignificantTrial = tfStimZScoresTrials > 6;
   mbReliableResponses = sum(tbSignificantTrial, 3) >= (size(tbSignificantTrial, 3)/2);
   %vbSignificantResponse = nanmax(mfStimZScores, [], 2) > 6;
   vbSignificantResponse = mfStimMeanResponses > 0.25;
   %vbResponsiveROI = vbSignificantResponse & any(mbReliableResponses, 2);
   vbResponsiveROI = sum(vbSignificantResponse.*mbReliableResponses,2) > 0;

   % - Save results
   save(strAnalysisFile, 'sRegions', 'mfStimZScores', 'vfBlankStds', ...
      'mfStimMeanResponses', 'mfStimStds', 'mfRegionTraces', 'tnTrialSampleSizes', ...
      'tfTrialResponses', 'cfTrialTraces', 'sNeuropilZThresh', 'tbSignificantTrial', 'mbReliableResponses', ...
      'vbSignificantResponse', 'vbResponsiveROI', '-append');
end


%% -- Plot figures

if (bPlotFigures)
   [~, vnRespROIs] = sort(nanmax(mfStimZScores, [], 2), 1, 'descend');
   if (nnz(vbResponsiveROI) < 10)
      vnRespROIs = vnRespROIs(1:10);
   else
      vnRespROIs = vnRespROIs(vbResponsiveROI(vnRespROIs));
   end   
   
   vhFigures = figure;
   MakeROIOverviewFigure(fsStack, fsStack.fPixelsPerUM, 1:2, sRegions, vnRespROIs, true);
   
   % - Make a plot of the response traces for different stimuli
   if (numel(vnRespROIs) > 0)
      vnPlotROIs = vnRespROIs(1:min(numel(vnRespROIs), 40));
      if (numel(vnRespROIs) == 1)
         vnPlotROIs = vnPlotROIs * [1 1];
      end
      hFigure2 = PlotRegionResponses(mfRegionTraces, fsStack, sRegions, mtStimLabelTimes, 0.1, 5, vnPlotROIs, {}, [], [], tBaseTimeShift);
      vhFigures = [vhFigures hFigure2];
   end
else
   vhFigures = [];
end


%% -- Store analysis data

% - Make a mask matrix
mbResponseMask = (labelmatrix(sRegions) > 0);

sAnalysis.sCellRegions = sRegions;
sAnalysis.mbResponseMask = mbResponseMask;
sAnalysis.vfBlankStds = vfBlankStds;
sAnalysis.mfStimMeans = mfStimMeanResponses;
sAnalysis.mfStimStds = mfStimStds;
sAnalysis.mfStimZScores = mfStimZScores;
sAnalysis.tnTrialSampleSizes = tnTrialSampleSizes;
sAnalysis.sNeuropilZThresh = sNeuropilZThresh;
sAnalysis.tfTrialResponses = tfTrialResponses;
sAnalysis.tbSignificantTrial = tbSignificantTrial;
sAnalysis.mbReliableResponses = mbReliableResponses;
sAnalysis.vbSignificantResponse = vbSignificantResponse;
sAnalysis.vbResponsiveROI = vbResponsiveROI;

end

%% AssignBaselineDistribution - FUNCTION Assign the baseline response for each stimulus

function fsStack = AssignBaselineDistribution(fsStack, vnUseStimIDs, nBlankStimID, tBaseTimeShift, fBaselineProportion, mtBlankTimes)

% -- Assign baseline
disp('--- AnalyseExpSetPeaks/AssignBaselineDistribution: Assigning baseline...');

% - Get an alignment mask
mbAlignMask = fsStack.GetAlignedMask; %#ok<NASGU>

% - Do not use baseline shift for assigning blank period, to avoid
% contamination with rising transients
[vtGlobalTime, ...
   vnBlockIndex, ~, ~, ...
   vnStimulusSeqID, vtTimeInStimPresentation, ...
   ~, vbUseFrame] = ...
   FrameStimulusInfo(fsStack, 1:size(fsStack, 3), 0);

[~, vnBlockIndexShift, ~, ~, ...
   vnStimulusSeqIDShift] = ...
   FrameStimulusInfo(fsStack, 1:size(fsStack, 3), tBaseTimeShift);


nProgress = 0;
nProgressMax = numel(fsStack.cstrFilenames) * numel(vnUseStimIDs);

fprintf(1, 'Assigning baseline: %6.2f%%', 0);

for (nBlock = 1:numel(fsStack.cstrFilenames))
   vbBlockFrames = vnBlockIndex == nBlock;
   vbBlockFramesShift = vnBlockIndexShift == nBlock;
   
   % - Is there a blank stimulus?
   if (~isempty(nBlankStimID))
      % - Get the baseline for this block from the blank stimulus
      vbBlockBlankStimFrames = vbBlockFrames & (vnStimulusSeqID == nBlankStimID) & vbUseFrame;
      
      % - Extract baseline frames
      tfBlankFrames = sort(double(fsStack.AlignedStack(:, :, vbBlockBlankStimFrames, 1)), 3, 'ascend');
      
      % - Find lower threshold
      nThrehsoldInd = floor(fBaselineProportion * nnz(vbBlockBlankStimFrames));
      mfThisBaseline = tfBlankFrames(:, :, nThrehsoldInd);
      
     
      
      % - Filter zeros, replace with NaN
      mfThisBaseline(mfThisBaseline <= 0) = nan;
      
      % - Calculate std. dev. for divisive normalisation
      mfBlockBlankStd = nanstd(bsxfun(@rdivide, tfBlankFrames, mfThisBaseline), [], 3);
      
      % - Assign the blank frame to the whole block (by default)
      fsStack.AssignBlankFrame(cat(3, mfThisBaseline, mfBlockBlankStd), vbBlockFrames);
      
   else
      % - Use all pre-stim blank periods
      vbAllPresBlankFrames = false(size(vtGlobalTime));
      for (nStimSeqID = vnUseStimIDs)
         vbAllPresBlankFrames =  vbAllPresBlankFrames | (vbBlockFrames & ...
                                    (vnStimulusSeqID == nStimSeqID) & vtTimeInStimPresentation >= mtBlankTimes(nStimSeqID, 1) & ...
                                    (vnStimulusSeqID == nStimSeqID) & vtTimeInStimPresentation <= mtBlankTimes(nStimSeqID, 2));
      end
      
      % - Extract baseline frames
      tfBlankFrames = sort(double(fsStack.AlignedStack(:, :, vbAllPresBlankFrames, 1)), 3, 'ascend');
      
      % - Find lower threshold
      nThrehsoldInd = floor(fBaselineProportion * nnz(vbAllPresBlankFrames));
      mfThisBaseline = tfBlankFrames(:, :, nThrehsoldInd);

      % - Filter zeros, replace with NaN
      mfThisBaseline(mfThisBaseline <= 0) = nan;
      
      % - Calculate std. dev. for divisive normalisation
      mfBlockBlankStd = nanstd(bsxfun(@rdivide, tfBlankFrames, mfThisBaseline), [], 3);
      
      % - Assign the blank frame to the whole block (by default)
      fsStack.AssignBlankFrame(cat(3, mfThisBaseline, mfBlockBlankStd), vbBlockFrames);
   end
   
   % - Assign separate blank frames for each other stimulus segment
   for (nStimSeqID = vnUseStimIDs)
      % - Find all frames for this presentation
      vbPresFrames = vbBlockFramesShift & (vnStimulusSeqIDShift == nStimSeqID);
      vbBlankFrames = vbBlockFrames & (vnStimulusSeqID == nStimSeqID) & ...
                      vtTimeInStimPresentation >= mtBlankTimes(nStimSeqID, 1) & ...
                      vtTimeInStimPresentation <= mtBlankTimes(nStimSeqID, 2);
     
      % - Estimate lower threshold
      tfPresBlank = sort(double(fsStack.AlignedStack(:, :, vbBlankFrames, 1)), 3, 'ascend');
      nThrehsoldInd = floor(fBaselineProportion * nnz(vbBlankFrames));
      mfThisBaseline = tfPresBlank(:, :, nThrehsoldInd);

      % - Filter zeros, replace with NaN
      mfThisBaseline(mfThisBaseline <= 0) = nan;

      % - Assign only the mean; Std.Dev is taken from block blank
      fsStack.AssignBlankFrame(cat(3, mfThisBaseline, mfBlockBlankStd), vbPresFrames);
      
      nProgress = nProgress + 1;
      fprintf(1, '\b\b\b\b\b\b\b%6.2f%%', nProgress / nProgressMax * 100);
   end
end

fprintf(1, '\n');

end

