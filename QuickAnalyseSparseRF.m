function QuickAnalyseSparseRF(cstrFilenames, ...
   vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
   tPixelDuration, tPixelBlankTime, tBlankStimTime, ...
   bAlign, bAssignBlack, bAssignBlank, ...
   tForceStackFrameDuration, cvnForceStackSequenceIDs, strImageJRoiSet, ...
   mfCustomAlignment, bOrderROIs, fhExtractionFunction)



% QuickAnalyseSparseRF - FUNCTION
%
% Usage: QuickAnalyseSparseRF(  cstrFilenames, ...
%                               vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
%                               tPixelDuration, tPixelBlankTime, tBlankStimTime, ...
%                               bAlign, bAssignBlack, bAssignBlank, ...
%                               tForceStackFrameDuration, cvnForceStackSequenceIDs, strImageJRoiSet, ...
%                               mfCustomAlignment, bOrderROIs, fhExtractionFunction)

% -- Defaults

DEF_bAlign = false;
DEF_bAssignBlack = false;
DEF_bAssignBlank = false;
DEF_mfCustomAlignment = 10;
DEF_bOrderROIs = true;
DEF_fhExtractionFunction = ExtractMean(1, bAssignBlank);

if (~exist('fhExtractionFunction', 'var') || isempty(fhExtractionFunction))
   fhExtractionFunction = DEF_fhExtractionFunction;
end

%% -- Check arguments, assign defaules

if (nargin < 8)
   disp('*** QuickAnalyseSparseRF: Incorrect usage');
   help QuickAnalyseSparseRF;
   return;
end

if (~exist('bAlign', 'var') || isempty(bAlign))
   bAlign = DEF_bAlign;
end

if (~exist('bAssignBlack', 'var') || isempty(bAssignBlack))
   bAssignBlack = DEF_bAssignBlack;
end

if (~exist('bAssignBlank', 'var') || isempty(bAssignBlank))
   bAssignBlank = DEF_bAssignBlank;
end

if (~exist('mfCustomAlignment', 'var') || isempty(mfCustomAlignment))
   mfCustomAlignment = DEF_mfCustomAlignment;
end

if (~exist('bOrderROIs', 'var') || isempty(bOrderROIs))
   bOrderROIs = DEF_bOrderROIs;
end


%% -- Generate a consistent "temporary" filename

sConfig.cstrFilenames = cstrFilenames;
sConfig.bAlign = bAlign;
sConfig.bAssignBlack = bAssignBlack;
sConfig.bAssignBlank = bAssignBlank;

strTempFilename = fullfile(tempdir, sprintf('QASRF_analysis_%u.mat', HashStructure(sConfig)));


%% -- Try to make a focus stack

try
   % - Construct a stack
   disp('--- QuickAnalyseSparseRF: Creating FocusStack...');
   fsStack = FocusStack(cstrFilenames);
   
   % - Assign frame rate, if required
   if (exist('tForceStackFrameDuration', 'var') && ~isempty(tForceStackFrameDuration))
      fsStack.tFrameDuration = tForceStackFrameDuration;
   end
   
   % - Assign stimulus sequence, if required
   if (exist('cvnForceStackSequenceIDs', 'var') && ~isempty(cvnForceStackSequenceIDs))
      if (~iscell(cvnForceStackSequenceIDs))
         cvnForceStackSequenceIDs = repmat({cvnForceStackSequenceIDs}, numel(cstrFilenames), 1);
      end
      
      fsStack.cvnSequenceIDs = cvnForceStackSequenceIDs;
   end
   
   % - Assign stimulus durations
   vtStimDurations = repmat(tPixelBlankTime + tPixelDuration, prod(vnNumPixels), 1);
   
   if (~isempty(tBlankStimTime) && (tBlankStimTime ~= 0))
      vtStimDurations = [tBlankStimTime; vtStimDurations];
   end
   
   fsStack.vtStimulusDurations = vtStimDurations;
   
   % - Assign stimulus use times
   mtStimUseTimes = repmat([tPixelBlankTime tPixelBlankTime+tPixelDuration], prod(vnNumPixels), 1);
   
   if (~isempty(tBlankStimTime) && (tBlankStimTime ~= 0))
      mtStimUseTimes = [[0 tBlankStimTime]; mtStimUseTimes];
   end
   
   fsStack.mtStimulusUseTimes = mtStimUseTimes;
   
   % -- Align stack, if requested
   if (bAlign)
      disp('--- QuickAnalyseSparseRF: Aligning...');
      nAlignChannel = size(fsStack, 4);
      fsStack.Align(nAlignChannel, false, 10, mfCustomAlignment);
   end
   
   % -- Assign black, if requested
   if (bAssignBlack)
      disp('--- QuickAnalyseSparseRF: Assigning black region...');
      fsStack.DefineBlackRegion;
   end
   
   % -- Assign blank frames, if requested
   if (bAssignBlank)
      disp('--- QuickAnalyseSparseRF: Assigning blank frames...');
      AssignBlank(fsStack, fhExtractionFunction);
   end
   
catch meErr
   % - Try to save the stack, if possible
   if (exist('fsStack', 'var'))
      save(strTempFilename, 'fsStack', 'vnNumPixels', 'fPixelOverlap', 'fPixelSizeDeg', 'vfScreenSizeDeg', 'tBlankStimTime');
   end
   
   rethrow(meErr);
end

fsStack.fPixelsPerUM = 2*fsStack.fPixelsPerUM;
% fsStack.fPixelsPerUM = 5;

% - Save stack
save(strTempFilename, 'fsStack', 'vnNumPixels', 'fPixelOverlap', 'fPixelSizeDeg', 'vfScreenSizeDeg', 'tBlankStimTime');


%% -- Define ROIs

% -- Automatically define regions
if (~exist('strImageJRoiSet', 'var') || isempty(strImageJRoiSet))   
   disp('--- QuickAnalyseSparseRF: Identifying ROIs...');
   if (size(fsStack, 4) > 1)
      sRegions = FindCells_GRChannels(fsStack);
   else
      sRegions = FindCells_GChannel(fsStack);
   end
else
   sRegions = ROIs2Regions(ReadImageJROI(strImageJRoiSet), size(fsStack, [1 2]));
end

% - Add a "neuropil" region
nNumPixels = prod(size(fsStack, 1:2)); %#ok<PSIZE>
sRegionsPlusNP = sRegions;
% sRegionsPlusNP.NumObjects = sRegionsPlusNP.NumObjects+1;
% sRegionsPlusNP.PixelIdxList = [{1:nNumPixels} sRegionsPlusNP.PixelIdxList];

% - Save ROIs
save(strTempFilename, 'sRegions', '-append');


%% -- Extract responses, order by Z-score

% - Extract responses, order by Z-score
disp('--- QuickAnalyseSparseRF: Extracting responses...');

[vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, tfTrialResponses, tnTrialSampleSizes] = ...
   ExtractRegionResponses(fsStack, sRegionsPlusNP, 1, fhExtractionFunction);

if (bAssignBlank)
   disp('--- QuickAnalyseSparseRF: Measuring responsiveness...');
   % - Compute Z-score against blank STDs
   nNumStimuli = numel(fsStack.vtStimulusDurations);
   mnStimSampleSizes = sum(tnTrialSampleSizes, 3);
   mfBlankStdsCorr = repmat(vfBlankStds', 1, nNumStimuli) ./ sqrt(mnStimSampleSizes);
   mfStimZScores = mfStimMeanResponses ./ mfBlankStdsCorr;
   
   tfBlankStdsCorr = repmat(vfBlankStds', [1, nNumStimuli, size(tfTrialResponses, 3)]) ./ sqrt(tnTrialSampleSizes);
   tfStimZScoresTrials = tfTrialResponses ./ tfBlankStdsCorr;
   
   % - Sort by Z-score
   if (bOrderROIs)
      vfMaxZScores = max(mfStimZScores, [], 2);   
      vfMaxZScores(1) = inf;  % - Force NP inclusion, as first region
      [nul, vnSort] = sort(vfMaxZScores, 'descend');
      
      sRegionsPlusNP.PixelIdxList = sRegionsPlusNP.PixelIdxList(vnSort);
      sRegionsPlusNP.vfRegionMaxZScores = vfMaxZScores(vnSort);
      mfRegionTraces = mfRegionTraces(vnSort, :);
      tfTrialResponses = tfTrialResponses(vnSort, :, :); %#ok<NASGU>
   end
   
else
   mfStimZScores = [];
   tfBlankStdsCorr = [];
   tfStimZScoresTrials = [];
end

% - Save responses
save(strTempFilename, 'vfBlankStds', 'mfStimMeanResponses', 'mfStimStds', ...
   'mfRegionTraces', 'tfTrialResponses', 'tnTrialSampleSizes', 'sRegionsPlusNP', ...
   'mfStimZScores', 'tfBlankStdsCorr', 'tfStimZScoresTrials', '-append');


%% -- Create figures

% % - Create overview figure
% nNumROIs = numel(sRegionsPlusNP.PixelIdxList);
% vhFigures = ...
%    MakeROIOverviewFigure(fsStack, [], 1:size(fsStack, 4), sRegionsPlusNP, 2:nNumROIs, true);
% 
% 
% % - Generate a mesh
% vfScreenX = 0:(1/DEF_nPointsPerDeg):vfScreenSizeDeg(1);
% vfScreenY = 0:(1/DEF_nPointsPerDeg):vfScreenSizeDeg(2);

RF_explorer(fsStack, vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
      vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, ...
      tfTrialResponses, tnTrialSampleSizes, sRegionsPlusNP, tBlankStimTime, ...
      mfStimZScores, tfBlankStdsCorr, tfStimZScoresTrials);



% --- END of QuickAnalyseSparseRF FUNCTION ---

%% AssignBlank FUNCTION - Assign blank frames to the stack
%
% Usage: AssignBlank(fsStack, fhExtractionFunction)

   function AssignBlank(fsStack, fhExtractionFunction)
      
      mbAlignMask = fsStack.GetAlignedMask; %#ok<NASGU>
      
      [vtGlobalTime, ...
         vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
         vnStimulusSeqID, vtTimeInStimPresentation, ...
         nul, vbUseFrame] = ...
         FrameStimulusInfo(fsStack, 1:size(fsStack, 3));
      
      % - Turn off blank warning and unaligned warning
      wOld = warning('off', 'FocusStack:FilteringBlankFrames');
      warning('off', 'FocusStack:BlankContainsNaN');
      warning('off', 'FocusStack:UnalignedStack');
      
      % - Turn off DFF conversion, just in case
      fsStack.BlankNormalisation('none');
      
      for (nBlock = 1:numel(fsStack.cstrFilenames)) %#ok<FORPF>
         % - Get the median blank frame for this block
         vbBlockFrames = vnBlockIndex == nBlock;
         vbBlockBlankStimFrames = vbBlockFrames & (vnStimulusSeqID == 1);
         vbBlockBlankFrames = vbBlockBlankStimFrames & vbUseFrame;
         
         tfBlankFrames = reshape(fhExtractionFunction(fsStack, :, vbBlockBlankFrames), size(fsStack, 1), size(fsStack, 2), []);
         mfThisBlankAvg = nanmedian(tfBlankFrames, 3);
         
         % - Filter zeros, replace with NaN
         mfThisBlankAvg(mfThisBlankAvg == 0) = nan;
         
         % - Calculate std. dev. for divisive normalisation
         mfThisBlankStd = nanstd(tfBlankFrames ./ repmat(mfThisBlankAvg, [1 1 nnz(vbBlockBlankFrames)]), [], 3);
         
         % - Assign the blank frame to the whole block (by default)
         fsStack.AssignBlankFrame(cat(3, mfThisBlankAvg, mfThisBlankStd), vbBlockFrames);
         
         % - Assign separate blank frames for each other stimulus segment
         for (nStimSeqID = unique(vnStimulusSeqID))
            % - Find all frames for this presentation
            vbPresFrames = vbBlockFrames & (vnStimulusSeqID == nStimSeqID);
            
            % - Find blank frames for this presentation (start and end of blank time)
            vbPresBlankFrames =  vbPresFrames & ...
               (vtTimeInStimPresentation >= 0) & ...
               (vtTimeInStimPresentation <= tPixelBlankTime);
            
            % - Extract these blank frames
            tfPresBlank = reshape(fhExtractionFunction(fsStack, :, vbPresBlankFrames), size(fsStack, 1), size(fsStack, 2), []);
            
            % - Assign the mean; Std.Dev is taken from block blank
            fsStack.AssignBlankFrame(cat(3, nanmedian(tfPresBlank, 3), mfThisBlankStd), vbPresFrames);
         end
      end
      
      % - Restore warnings
      warning(wOld);
      
   end


end
% --- END of QUickAnalyseSparseRF.m ---
