function [sRFAnalysis] = OfflineAnalysisSparseRF(cstrFilenames, ...
   vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
   tPixelDuration, nBlankStimID, tBlankStimTime, ...
   bAlign, bAssignBlack, bAssignBlank, ...
   tForceStackFrameDuration, cvnForceStackSequenceIDs, strImageJRoiSet, ...
   mfCustomAlignment, fhExtractionFunction)



% OfflineAnalysisSparseRF - FUNCTION
%
% Usage: [sRFAnalysis] = OfflineAnalysisSparseRF(  cstrFilenames, ...
%                               vnNumPixels, fPixelOverlap, fPixelSizeDeg, vfScreenSizeDeg, ...
%                               tPixelDuration, nBlankStimID, tBlankStimTime, ...
%                               bAlign, bAssignBlack, bAssignBlank, ...
%                               tForceStackFrameDuration, cvnForceStackSequenceIDs, strImageJRoiSet, ...
%                               mfCustomAlignment, fhExtractionFunctionFunction)
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
   disp('*** OfflineAnalysisSparseRF: Incorrect usage');
   help OfflineAnalysisSparseRF;
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

if (~isa(cstrFilenames, 'FocusStack'))
   % - Construct a stack
   disp('--- QuickAnalyseSparseRF: Creating FocusStack...');
   fsStack = FocusStack(cstrFilenames);
else
   fsStack = cstrFilenames;
end

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

% fsStack.fPixelsPerUM = 2*fsStack.fPixelsPerUM;
% fsStack.fPixelsPerUM = 5;

% - Save stack
sRFAnalysis.fsStack = fsStack;
sRFAnalysis.vnNumPixels = vnNumPixels;
sRFAnalysis.fPixelOverlap = fPixelOverlap;
sRFAnalysis.fPixelSizeDeg = fPixelSizeDeg;
sRFAnalysis.vfScreenSizeDeg = vfScreenSizeDeg;
sRFAnalysis.tBlankStimTime = tBlankStimTime;


%% -- Define ROIs

% -- Automatically define regions
if (~exist('strImageJRoiSet', 'var') || isempty(strImageJRoiSet))   
   disp('--- QuickAnalyseSparseRF: Identifying ROIs...');
   if (size(fsStack, 4) > 1)
      sRegions = FindCells_GRChannels(fsStack);
   else
      sRegions = FindCells_GChannel(fsStack);
   end
   
elseif (ischar(strImageJRoiSet))
   sRegions = ROIs2Regions(ReadImageJROI(strImageJRoiSet), size(fsStack, [1 2]));
   
else
   sRegions = strImageJRoiSet;
end

sRFAnalysis.sRegions = sRegions;


%% -- Extract responses, order by Z-score

% - Extract responses, order by Z-score
disp('--- QuickAnalyseSparseRF: Extracting responses...');

[vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, tfTrialResponses, tnTrialSampleSizes] = ...
   ExtractRegionResponses(fsStack, sRFAnalysis.sRegions, nBlankStimID, fhExtractionFunction(bAssignBlank));

if (bAssignBlank)
   disp('--- QuickAnalyseSparseRF: Measuring responsiveness...');
   % - Compute Z-score against blank STDs
   mnStimSampleSizes = sum(tnTrialSampleSizes, 3);
   mfBlankStdsCorr = repmat(vfBlankStds', 1, nNumStimuli) ./ sqrt(mnStimSampleSizes);
   mfStimZScores = mfStimMeanResponses ./ mfBlankStdsCorr;
   
   tfBlankStdsCorr = repmat(vfBlankStds', [1, nNumStimuli, size(tfTrialResponses, 3)]) ./ sqrt(tnTrialSampleSizes);
   tfStimZScoresTrials = tfTrialResponses ./ tfBlankStdsCorr;
      
else
   mfStimZScores = [];
   tfBlankStdsCorr = [];
   tfStimZScoresTrials = [];
end

sRFAnalysis.vfBlankStds = vfBlankStds;
sRFAnalysis.mfStimMeanResponses = mfStimMeanResponses;
sRFAnalysis.mfStimStds = mfStimStds;
sRFAnalysis.mfRegionTraces = mfRegionTraces;
sRFAnalysis.tfTrialResponses = tfTrialResponses;
sRFAnalysis.tnTrialSampleSizes = tnTrialSampleSizes;
sRFAnalysis.mfStimZScores = mfStimZScores;
sRFAnalysis.tfBlankStdsCorr = tfBlankStdsCorr;
sRFAnalysis.tfStimZScoresTrials = tfStimZScoresTrials;


%% -- Extract all RFs

% -- Make RF export structures

fprintf(1, '--- OfflineAnalysisSparseRF: Exporting RFs for all ROIs [%6.2f%%]', 0);

handles = [];

for (nROIIndex = sRFAnalysis.sRegions.NumObjects:-1:1)
   % - Get this RF export structure
   [vsRFExport(nROIIndex), handles] = RFExportStruct(handles, sRFAnalysis, nROIIndex);
   
   fprintf(1, '\b\b\b\b\b\b\b\b%6.2f%%]', (1 - nROIIndex / sRFAnalysis.sRegions.NumObjects) * 100);
end

fprintf(1, '\b\b\b\b\b\b\b\b%6.2f%%]\n', 100);

sRFAnalysis.vsRFEstimates = vsRFExport;


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

%% MakeRFImage FUNCTION -- Estimate RFs for a set of ROIs

   function [mfRFImage, handles] = MakeRFImage(handles, sRFAnalysis, vnSelectedROIs)
      % - Defaults
      DEF_nRFSamplesPerDeg = 4;      
      
      % - Make a mesh
      if (~isfield(handles, 'mfXMesh'))
         vfX = 0:(1/DEF_nRFSamplesPerDeg):sRFAnalysis.vfScreenSizeDeg(1);
         vfY = 0:(1/DEF_nRFSamplesPerDeg):sRFAnalysis.vfScreenSizeDeg(2);
         [handles.mfXMesh, handles.mfYMesh] = ndgrid(vfX, vfY);
      end
      
      % - Work out stimulus locations
      if (~isfield(handles, 'tfGaussian'))
         vfStimSizeDeg = sRFAnalysis.fPixelSizeDeg .* sRFAnalysis.vnNumPixels * (1-sRFAnalysis.fPixelOverlap);
         vfXCentres = linspace(-vfStimSizeDeg(1)/2, vfStimSizeDeg(1)/2, sRFAnalysis.vnNumPixels(1) + 1) + (sRFAnalysis.fPixelSizeDeg * (1-sRFAnalysis.fPixelOverlap))/2;
         vfXCentres = vfXCentres(1:end-1);
         vfYCentres = linspace(-vfStimSizeDeg(2)/2, vfStimSizeDeg(2)/2, sRFAnalysis.vnNumPixels(2) + 1) + (sRFAnalysis.fPixelSizeDeg * (1-sRFAnalysis.fPixelOverlap))/2;
         vfYCentres = vfYCentres(1:end-1);
         
         vfXStimCentres = vfXCentres + sRFAnalysis.vfScreenSizeDeg(1)/2;
         vfYStimCentres = vfYCentres + sRFAnalysis.vfScreenSizeDeg(2)/2;
         
         [handles.mfXStimCentreMesh, handles.mfYStimCentreMesh] = meshgrid(vfXStimCentres, vfYStimCentres);
         
         % - Pre-compute distance meshes and unitary Gaussians
         for (nStimID = prod(sRFAnalysis.vnNumPixels):-1:1)
            tfDistanceMeshSqr(:, :, nStimID) = (handles.mfXMesh - handles.mfXStimCentreMesh(nStimID)).^2 + (handles.mfYMesh - handles.mfYStimCentreMesh(nStimID)).^2;
         end
         
         fRFSigma = sRFAnalysis.fPixelSizeDeg/4 * 2;
         handles.tfGaussian = exp(-1/(2*fRFSigma.^2) .* tfDistanceMeshSqr);
      end
      
      % - Iterate over ROIs and build up a Gaussian RF estimate
      mfRFImage = zeros(size(handles.mfXMesh));
      
      % - Should we use Z-scored responses, or raw responses?
      if (~isempty(sRFAnalysis.mfStimZScores))
         bUsedFF = true;
      else
         bUsedFF = false;
      end
      
      % - Get the responses
      if (bUsedFF)
         mfStimResp = sRFAnalysis.mfStimZScores(vnUseStimIDs, :);
         fThreshold = 3;
      else
         mfStimResp = sRFAnalysis.mfStimMeanResponses(vnUseStimIDs, :);
         fThreshold = 0;
      end
      
      % - Clip below threshold responses
      mfStimResp(mfStimResp < fThreshold) = nan;
      
      for (MRFI_nROIIndex = 1:numel(vnSelectedROIs))
         % - Accumulate Gaussians over stimulus locations for this RF
         nROI = vnSelectedROIs(MRFI_nROIIndex);
         
         vfROIResponse = mfStimResp(nROI, :);
         vfROIResponse = permute(vfROIResponse, [3 1 2]);
         
         % - Remove NaNs and INFs
         vfROIResponse(isnan(vfROIResponse)) = 0;
         vfROIResponse(isinf(vfROIResponse)) = 0;
         
         % - Estimate RF
         mfRFImage = mfRFImage + sum(bsxfun(@times, handles.tfGaussian, vfROIResponse), 3) ./ numel(vnSelectedROIs);         
      end
      
      if (bUsedFF)
         handles.strRFDataDescription = 'RF values are the Z-score of response';
      else
         handles.strRFDataDescription = 'RF values are the raw response';
      end
   end

%% RFExportStruct FUNCTION - Create an export struct for an RF estimate

   function [sExport, handles] = RFExportStruct(handles, sRFAnalysis, vnSelectedROIs)
      
      % -- Create data for export
      sExport.vnSelectedROIs = vnSelectedROIs;
      [sExport.mfRFImage, handles] = MakeRFImage(handles, sRFAnalysis, sExport.vnSelectedROIs);
      
      % -- Find population RF center and size
      
      if (~isempty(sRFAnalysis.mfStimZScores))
         fRFThreshold = 3;
      else
         fRFThreshold = eps;
      end
      
      mnRFPeaks = FindPeaks(sExport.mfRFImage, ...
         sRFAnalysis.fPixelSizeDeg .* size(sExport.mfRFImage, 1) ./ sRFAnalysis.vfScreenSizeDeg(1), ...
         sExport.mfRFImage > fRFThreshold);
      
      % - Find maximum peak
      vfRFMagnitude = sExport.mfRFImage(sub2ind(size(sExport.mfRFImage), mnRFPeaks(:, 2), mnRFPeaks(:, 1)));
      [nul, nRFIndex] = nanmax(vfRFMagnitude);
      
      % - Convert to degrees offset from center of screen
      if (~isempty(mnRFPeaks))
         sExport.vfRFOffsetDeg = (mnRFPeaks(nRFIndex, [2 1]) - size(sExport.mfRFImage)./2) .* (sRFAnalysis.vfScreenSizeDeg ./ size(sExport.mfRFImage));
      else
         sExport.vfRFOffsetDeg = [];
      end
      
      % - Estimate population RF size
   end



end

% --- END of OfflineAnalysisSparseRF.m ---
