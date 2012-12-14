function [sNeuropilZThresholds] = ...
   ThresholdNeuropilResponse(fsData, sRegions, vnUseStimulusSeqIDs, fAlpha, nBlankStimID, fhExtractionFunction)

% ThresholdNeuropilResponse - FUNCTION Estimate a threshold for rejecting neuropil contaminated responses
%
% Usage: [sNeuropilZThresholds] = ...
%    ThresholdNeuropilResponse(fsData, sRegions, vnUseStimulusSeqIDs, fAlpha, nBlankStimID, fhExtractionFunction)
% Usage: [sNeuropilZThresholds] = ...
%    ThresholdNeuropilResponse(fsData, sRegions, vnUseStimulusSeqIDs, fAlpha, vbBlankFrames, fhExtractionFunction)
%
% 'fsData' is a FocusStack object, with blank frames assigned.
%
% 'sRegions' is a regions structure, as returned by bwconncomp, that
% defines the ROIs of identified cells.  "Neuropil" is defined by the area
% OUTSIDE these regions.
%
% 'vnUseStimulusSeqIDs' defines the stimulus sequence IDs that correspond
% to experimental conditions with an activity-evoking stimulus, as opposed
% to "blank" stimuli.
%
% 'fAlpha' defines the alpha threshold above which responses will be deemed
% significantly above neuropil.
%
% The optional argument 'fhExtractionFunction' is a function handle that
% extracts responses from the stack.  It must be a function handle with the
% signature [mfRawTrace, varagout] = fh(fsData, vnPixels, vnFrames),
% identical to that required by ExtractRegionResponses.
%
% 'sNeuropilZThresholds' will be a matlab structure with the folowing
% fields:
%    .fSingleStimMean - The Z-score threshold that applies to a single
%       stimulus, averaged over the duration of a stimulus presentation
%       period, FOR A SINGLE PIXEL.  This value must be corrected by
%       multiplying by sqrt(nNumPixels) to compare with the average
%       response of a whole ROI.
%    .fSingleStimResp - The Z-score threshold that applies to a single
%       stimulus, if the "response" is considered as the peak respose
%       during a stimulus presentation.  This is a per-pixel threshold,
%       normliased for the number of samples reported by the extraction
%       function.
%    .fMaxStimMean - The Z-score threshold that applies to the max Z-score
%       over all stimuli defined in 'fsData', when each "response" is the
%       average over the duration of a stimulus.  This is a per-pixel
%       threshold, and must be corrected to compare with a ROI average.
%    .fMaxStimResp - The Z-score threshold that applies to the max Z-score
%       over all stimuli, when each "response" is the extracted response
%       according to the extraction function.
%    .fAlpha - The alpha threshold used to derive these thresholds
%

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 24th October, 2011

% -- Check arguments

if (nargin < 4)
   disp('*** ThresholdNeuropilResponse: Incorrect usage');
   help ThresholdNeuropilResponse;
   return;
end

if (~exist('fhExtractionFunction', 'var') || isempty(fhExtractionFunction))
   fhExtractionFunction = ExtractMean;
end


% -- Get a list of matching stack frames

vnStackSize = size(fsData);

[vtGlobalTime, ...
 vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
 vnStimulusSeqID, vtTimeInStimPresentation, ...
 vnPresentationIndex, vbUseFrame] = ...
 	FrameStimulusInfo(fsData, 1:vnStackSize(3));

vbMatchingStackFrames = ismember(vnStimulusSeqID, vnUseStimulusSeqIDs);


% -- Check whether blank frames have been supplied

if (isempty(fsData.mnAssignedBlankMeanFrames) || isempty(fsData.mnAssignedBlankStdFrames) || ...
    any(isnan(fsData.mnAssignedBlankMeanFrames(vbMatchingStackFrames, 1))) || ...
    any(isnan(fsData.mnAssignedBlankStdFrames(vbMatchingStackFrames, 1))))
   error('FocusStack:BlankFramesRequired', ...
      '*** FocusStack/ResponsiveMask: Blank mean and std. dev. frames must be assigned to use ThresholdNeuropilResponse.');
end


% -- Find stimulus responsive pixels

% - Turn off blank normalisation
strOldNorm = fsData.BlankNormalisation('none');

% - Read blank data for neuropil
mbRegions = labelmatrix(sRegions);
mbNeuropil = ~mbRegions;
sNeuropilRegions = sRegions;
sNeuropilRegions.NumObjects = 1;
sNeuropilRegions.PixelIdxList = {find(mbNeuropil)};
nNumNeuropilPixels = nnz(mbNeuropil);

nNumStim = numel(vnUseStimulusSeqIDs);

% -- Extract blank frames
if (~isscalar(nBlankStimID))
   vbBlankFrames = nBlankStimID;
else
   vbBlankFrames = (vnStimulusSeqID == nBlankStimID) & vbUseFrame;
end

[mfNeuropilBlankTrace] = fhExtractionFunction(fsData, mbNeuropil(:), vbBlankFrames);
mfNeuropilBlankTrace = double(mfNeuropilBlankTrace);
vfNeuropilBlankMean = nanmean(mfNeuropilBlankTrace, 2);

vfBlankStds = nanstd(mfNeuropilBlankTrace, [], 2);

nNumBlocks = numel(fsData.cstrFilenames);

% - Find per-stim means and calculate Z score

for (nStim = nNumStim:-1:1)                  % Go backwards to pre-allocate
   % - Find matching frames, within desired stimulus use duration
   nThisStimID = vnUseStimulusSeqIDs(nStim);
   vbThisStimFrames = (vnStimulusSeqID == nThisStimID) & vbUseFrame;
   
   for (nBlock = nNumBlocks:-1:1)
      vbThisBlockStimFrames = vbThisStimFrames & (vnBlockIndex == nBlock);
      
      % - Extract response and compute mean
      [mfNeuropilStimTrace, nul, nul, nNumSamples, vfNPResponse] = fhExtractionFunction(fsData, mbNeuropil(:), vbThisBlockStimFrames);
      vfNeuropilMean = nanmean(double(mfNeuropilStimTrace), 2);

      % - Record responses and number of samples
      tfNPMeanResp(:, nStim, nBlock) = vfNeuropilMean;
      tfNPResp(:, nStim, nBlock) = vfNPResponse;
      
      mnMeanSamples(nStim, nBlock) = nnz(vbThisBlockStimFrames);
      mnRespSamples(nStim, nBlock) = nNumSamples;
   end
end

% - Compute corrected std.dev. for combined blank
% - Z score is (mean obs. - mean blank) / (std blank / sqrt(num stim frames))
%   This uses the expected standard error of the mean of several response
%   frames
mfNPMeanMeanStimZ = (nanmean(tfNPMeanResp, 3) - repmat(vfNeuropilBlankMean, 1, nNumStim)) ./ (repmat(vfBlankStds, 1, nNumStim) ./ sqrt(repmat(sum(mnMeanSamples, 2)', nNumNeuropilPixels, 1)));
mfNPRespMeanStimZ = (nanmean(tfNPResp, 3) - repmat(vfNeuropilBlankMean, 1, nNumStim)) ./ (repmat(vfBlankStds, 1, nNumStim) ./ sqrt(repmat(sum(mnRespSamples, 2)', nNumNeuropilPixels, 1)));

tfNPTrialMeanStimZ = (tfNPMeanResp - repmat(vfNeuropilBlankMean, [1, nNumStim, nNumBlocks])) ./ (repmat(vfBlankStds, [1, nNumStim, nNumBlocks]) ./ sqrt(repmat(permute(mnMeanSamples, [3 1 2]), [nNumNeuropilPixels, 1, 1])));
tfNPTrialRespStimZ = (tfNPResp - repmat(vfNeuropilBlankMean, [1, nNumStim, nNumBlocks])) ./ (repmat(vfBlankStds, [1, nNumStim, nNumBlocks]) ./ sqrt(repmat(permute(mnRespSamples, [3 1 2]), [nNumNeuropilPixels, 1, 1])));


% -- Determine significance thresholds for single stimulus responses

nNumSamples = numel(mfNPMeanMeanStimZ);
nThresholdIndex = ceil(nNumSamples * (1-fAlpha));

vfNPStimMeanZ = sort(mfNPMeanMeanStimZ(:));
vfNPStimRespZ = sort(mfNPRespMeanStimZ(:));

sNeuropilZThresholds.fAlpha = fAlpha;
sNeuropilZThresholds.fSingleStimMean = vfNPStimMeanZ(nThresholdIndex);
sNeuropilZThresholds.fSingleStimResp = vfNPStimRespZ(nThresholdIndex);


% -- Determine significance thresholds for the max of all stimuli

vfNeuropilMaxMeanStimZ = sort(nanmax(mfNPMeanMeanStimZ, [], 2));
vfNeuropilMaxRespStimZ = sort(nanmax(mfNPRespMeanStimZ, [], 2));

nNumSamples = numel(vfNeuropilMaxMeanStimZ);
nThresholdIndex = ceil(nNumSamples * (1-fAlpha));

sNeuropilZThresholds.fMaxStimMean = vfNeuropilMaxMeanStimZ(nThresholdIndex);
sNeuropilZThresholds.fMaxStimResp = vfNeuropilMaxRespStimZ(nThresholdIndex);


% -- Determine significance thresholds for single trials

vfNPStimTrialMeanZ = sort(tfNPTrialMeanStimZ(:));
vfNPStimTrialRespZ = sort(tfNPTrialRespStimZ(:));

nNumSamples = numel(vfNPStimTrialMeanZ);
nThresholdIndex = ceil(nNumSamples * (1-fAlpha));

sNeuropilZThresholds.fTrialStimMean = vfNPStimTrialMeanZ(nThresholdIndex);
sNeuropilZThresholds.fTrialStimResp = vfNPStimTrialRespZ(nThresholdIndex);


% - Restore blank normalisation 
fsData.BlankNormalisation(strOldNorm);

% --- END of ThresholdNeuropilResponse.m ---

