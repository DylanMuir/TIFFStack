function QuickAnalyseSparseRF(cstrFilenames, ...
   vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
   tPixelDuration, nBlankStimID, tBlankStimTime, ...
   bAlign, bAssignBlack, bAssignBlank, ...
   tForceStackFrameDuration, cvnForceStackSequenceIDs, strImageJRoiSet, ...
   mfCustomAlignment, bOrderROIs, fhExtractionFunction)



% QuickAnalyseSparseRF - FUNCTION
%
% Usage: QuickAnalyseSparseRF(  cstrFilenames, ...
%                               vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
%                               tPixelDuration, nBlankStimID, tBlankStimTime, ...
%                               bAlign, bAssignBlack, bAssignBlank, ...
%                               tForceStackFrameDuration, cvnForceStackSequenceIDs, strImageJRoiSet, ...
%                               mfCustomAlignment, bOrderROIs, fhExtractionFunctionFunction)
%
% 'mfCustomAlignment' can be a matrix, in which case it is used as an image
% against which to measure misalignment.  It can be a scalar, in which case
% is is used as a frame index to extract a reference image from the stack.
% Or is can be a function handle @(oStack)Align(oStack, params), which
% computes misalignment on a FocusStack with arbitrary parameters.

% -- Defaults

DEF_bAlign = false;
DEF_bAssignBlack = false;
DEF_bAssignBlank = false;
DEF_mfCustomAlignment = 10;
DEF_bOrderROIs = true;

fUsePreStimBlankProportion = 0.5;

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

DEF_fhExtractionFunction = @(bDFF)ExtractMean(1, bDFF);
if (~exist('fhExtractionFunction', 'var') || isempty(fhExtractionFunction))
   fhExtractionFunction = DEF_fhExtractionFunction;
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
   
   tPixelBlankTime = fsStack.tBlankTime;
   
   % - Assign stimulus durations and "use data" times
   nNumStimuli = prod(vnNumPixels) + ~isempty(nBlankStimID);
   
   vtStimulusDurations = repmat(tPixelDuration + tPixelBlankTime, nNumStimuli, 1);
   mtStimulusUseTimes = repmat(tPixelBlankTime + [0 tPixelDuration], nNumStimuli, 1);
   mtStimLabelTimes = repmat(tPixelBlankTime + [0 tPixelDuration], nNumStimuli, 1);
   mtBlankTimes = repmat(tPixelBlankTime * [fUsePreStimBlankProportion 1], nNumStimuli, 1);
   vnUseStimIDs = 1:nNumStimuli;
   
   if (~isempty(nBlankStimID))
       vtStimulusDurations(nBlankStimID) = tBlankStimTime;
       mtStimulusUseTimes(nBlankStimID, 1) = 0;
       mtStimLabelTimes(nBlankStimID, :) = nan;
       mtBlankTimes(nBlankStimID, :) = tBlankStimTime * [fUsePreStimBlankProportion 1];
       vnUseStimIDs = vnUseStimIDs(vnUseStimIDs ~= nBlankStimID);
   end
   
   fsStack.vtStimulusDurations = vtStimulusDurations;
   fsStack.mtStimulusUseTimes = mtStimulusUseTimes;
   
   % -- Align stack, if requested
   if (bAlign)
      disp('--- QuickAnalyseSparseRF: Aligning...');
      
      if (isa(mfCustomAlignment, 'function_handle'))
         fsStack.mfFrameShifts = mfCustomAlignment(fsStack);
      else
          
         nAlignChannel = size(fsStack, 4);
         fsStack.Align(nAlignChannel, false, 1, mfCustomAlignment);
      end
   end
   
   % -- Assign black, if requested
   if (bAssignBlack)
      disp('--- QuickAnalyseSparseRF: Assigning black region...');
      fsStack.DefineBlackRegion;
   end
   
   % -- Assign blank frames, if requested
   if (bAssignBlank)
      disp('--- QuickAnalyseSparseRF: Assigning blank frames...');
      AssignBlank(fsStack, fhExtractionFunction, mtBlankTimes, nBlankStimID);
   end
   
catch meErr
   % - Try to save the stack, if possible
   if (exist('fsStack', 'var'))
      save(strTempFilename, 'fsStack', 'vnNumPixels', 'fPixelOverlap', 'fPixelSizeDeg', 'vfScreenSizeDeg', 'tBlankStimTime');
      fprintf('*** QuickAnalyseSparseRF: Saved stack to [%s]\n', strTempFilename);
   end
   
   rethrow(meErr);
end

% fsStack.fPixelsPerUM = 2*fsStack.fPixelsPerUM;
% fsStack.fPixelsPerUM = 5

% - Save stack
save(strTempFilename, 'fsStack', 'vnNumPixels', 'fPixelOverlap', 'fPixelSizeDeg', 'vfScreenSizeDeg', 'nBlankStimID', 'tBlankStimTime');


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

% - Mask off mis-aligned regions
mbAlignedMask = fsStack.GetAlignedMask;
for (nRegion = 1:sRegionsPlusNP.NumObjects)
   sRegionsPlusNP.PixelIdxList{nRegion} = sRegionsPlusNP.PixelIdxList{nRegion}(mbAlignedMask(sRegionsPlusNP.PixelIdxList{nRegion}));
end

% - Save ROIs
save(strTempFilename, 'sRegions', '-append');


%% -- Extract responses, order by Z-score

% - Extract responses, order by Z-score
disp('--- QuickAnalyseSparseRF: Extracting responses...');

[vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, tfTrialResponses, tnTrialSampleSizes] = ...
   ExtractRegionResponses(fsStack, sRegionsPlusNP, nBlankStimID, fhExtractionFunction(bAssignBlank));

if (bAssignBlank)
   disp('--- QuickAnalyseSparseRF: Measuring responsiveness...');
   % - Compute Z-score against blank STDs
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
      tfTrialResponses, tnTrialSampleSizes, sRegionsPlusNP, nBlankStimID, tBlankStimTime, ...
      mfStimZScores, tfBlankStdsCorr, tfStimZScoresTrials);



% --- END of QuickAnalyseSparseRF FUNCTION ---

%% AssignBlank FUNCTION - Assign blank frames to the stack
%
% Usage: AssignBlank(fsStack, fhExtractionFunction)

   function AssignBlank(fsStack, fhExtractionFunction, mtBlankTimes, nBlankStimID)
      
      mbAlignMask = fsStack.GetAlignedMask; %#ok<NASGU>
      
      [vtGlobalTime, ...
         vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
         vnStimulusSeqID, vtTimeInStimPresentation, ...
         nul, vbUseFrame] = ...
         FrameStimulusInfo(fsStack, 1:size(fsStack, 3)); %#ok<SETNU>
      
      % - Turn off blank warning and unaligned warning
      wOld = warning('off', 'FocusStack:FilteringBlankFrames');
      warning('off', 'FocusStack:BlankContainsNaN');
      warning('off', 'FocusStack:UnalignedStack');
      
      % - Turn off DFF conversion, just in case
      fsStack.BlankNormalisation('none');
      fhEF = fhExtractionFunction(false);
      
      vnStimSet = unique(vnStimulusSeqID);
      vnStimSet = vnStimSet(~isnan(vnStimSet));
      vnStimSet = setdiff(vnStimSet, nBlankStimID);
      
      nProgress = 0;
      nProgressMax = numel(fsStack.cstrFilenames) * numel(vnStimSet);
      
      fprintf(1, 'Assigning blank: %6.2f%%', 0);
      
      for (nBlock = 1:numel(fsStack.cstrFilenames)) %#ok<FORPF>
         % - Get the median blank frame for this block
         vbBlockFrames = vnBlockIndex == nBlock;
         vbBlockBlankStimFrames = vbBlockFrames & (vnStimulusSeqID == nBlankStimID);
         vbBlockBlankFrames = vbBlockBlankStimFrames & vbUseFrame;
         
         tfBlankFrames = reshape(fhEF(fsStack, 1:prod(size(fsStack, 1:2)), vbBlockBlankFrames), size(fsStack, 1), size(fsStack, 2), []); %#ok<PSIZE>
         mfThisBlankAvg = nanmean(tfBlankFrames, 3);
         
         % - Filter zeros, replace with NaN
         mfThisBlankAvg(mfThisBlankAvg == 0) = nan;
         
         % - Calculate std. dev. for divisive normalisation
         mfThisBlankStd = nanstd(tfBlankFrames ./ repmat(mfThisBlankAvg, [1 1 nnz(vbBlockBlankFrames)]), [], 3);
         
         % - Assign the blank frame to the whole block (by default)
         fsStack.AssignBlankFrame(cat(3, mfThisBlankAvg, mfThisBlankStd), vbBlockFrames);
         
         % - Assign separate blank frames for each other stimulus segment
         for (nStimSeqID = vnStimSet)
            % - Find all frames for this presentation
            vbPresFrames = vbBlockFrames & (vnStimulusSeqID == nStimSeqID);
            
            % - Find blank frames for this presentation (start and end of blank time)
            vbPresBlankFrames =  vbPresFrames & ...
               (vtTimeInStimPresentation >= mtBlankTimes(nStimSeqID, 1)) & ...
               (vtTimeInStimPresentation <= mtBlankTimes(nStimSeqID, 2));
            
            % - Extract these blank frames
            tfPresBlank = reshape(fhEF(fsStack, 1:prod(size(fsStack, 1:2)), vbPresBlankFrames), size(fsStack, 1), size(fsStack, 2), []); %#ok<PSIZE>
            
            % - Assign the mean; Std.Dev is taken from block blank
            fsStack.AssignBlankFrame(cat(3, nanmean(tfPresBlank, 3), mfThisBlankStd), vbPresFrames);
            
            nProgress = nProgress + 1;
            fprintf(1, '\b\b\b\b\b\b\b%6.2f%%', nProgress / nProgressMax * 100);
         end
      end
      
      % - Restore warnings
      warning(wOld);
      
   end


end
% --- END of QUickAnalyseSparseRF.m ---
