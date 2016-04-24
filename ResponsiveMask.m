function [sRespRegions, mbRespMask, mfMaxStimZ, mfMaxStimZFilt, tfStimZFilt] = ...
   ResponsiveMask(oStack, vnUseStimulusSeqIDs, nChannel, fAlphaThreshold, nMinRegionSize, nMaxRegionSize)

% ResponsiveMask - METHOD Identify stimulus-responsive pixels in the stack
%
% Usage: [sRespRegions, mbRespMask, mfMaxZScores, mfMaxFiltZScores] = ...
%           ResponsiveMask(oStack, vnUseStimulusSeqIDs ...
%                          <, nChannel, fAlphaThreshold>, ...
%                          <nMinRegionSize, nMaxRegionSize>)
%
% 'oStack' is a FocusStack object.
%
% 'vnStimulusSeqIDs' is a list of stimulus sequence IDs to examine for
% responsivity. 'nChannel' is the signal channel to use (1 by default).
% 'fAlphaThreshold' is the threshold significance level (0.05 by default).

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 19th November, 2010

DEF_fAlphaThreshold = 0.05;

if (~exist('nChannel', 'var') || isempty(nChannel))
   nChannel = 1;
end

if (~exist('fAlphaThreshold', 'var') || isempty(fAlphaThreshold))
   fAlphaThreshold = DEF_fAlphaThreshold;
end

mtUseStimDurations = oStack.mtStimulusUseTimes;

if (~isscalar(nChannel))
   error('FocusStack:InvalidArgument', '*** FocusStack/ResponsiveMask: ''nChannel'' must be a scalar integer.');
end

if (isempty(vnUseStimulusSeqIDs))
   vnUseStimulusSeqIDs = 1:oStack.nNumStimuli;
end
   
vnUseStimulusSeqIDs = unique(vnUseStimulusSeqIDs(:));
if (any(vnUseStimulusSeqIDs > oStack.nNumStimuli) || any(vnUseStimulusSeqIDs < 1))
   error('FocusStack:InvalidArgument', ...
      '*** FocusStack/ResponsiveMask: ''vnStimulusSeqIDs'' must be limited to [1 .. %d] for this stack.', ...
      oStack.nNumStimuli);
end

vnStackSize = size(oStack);
if ((nChannel > vnStackSize(4)) || (nChannel < 1))
   error('FocusStack:InvalidArgument', ...
      '*** FocusStack/ResponsiveMask: ''nChannel'' must be limited to [1 .. %d] for this stack.', ...
      vnStackSize(4));
end


% -- Get a list of matching stack frames

[vtGlobalTime, ...
 vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
 vnStimulusSeqID, vtTimeInStimPresentation] = ...
 	FrameStimulusInfo(oStack, 1:vnStackSize(3));

vbMatchingStackFrames = ismember(vnStimulusSeqID, vnUseStimulusSeqIDs);


% -- Check whether blank frames have been supplied

if (isempty(oStack.mnAssignedBlankMeanFrames) || isempty(oStack.mnAssignedBlankStdFrames) || ...
    any(isnan(oStack.mnAssignedBlankMeanFrames(vbMatchingStackFrames, 1))) || ...
    any(isnan(oStack.mnAssignedBlankStdFrames(vbMatchingStackFrames, 1))))
   error('FocusStack:BlankFramesRequired', ...
      '*** FocusStack/ResponsiveMask: Blank mean and std. dev. frames must be assigned to use ResponsiveMask.');
end


% -- Turn off DFF conversion, if present

strNormalisation = oStack.BlankNormalisation('none');


% -- Find stimulus responsive pixels

nNumStim = numel(vnUseStimulusSeqIDs);

% - Find per-stim means and calculate Z score
for (nStim = nNumStim:-1:1)                  % Go backwards to pre-allocate
   % - Find matching frames, within desired stimulus use duration
   nThisStimID = vnUseStimulusSeqIDs(nStim);
   vbThisStimFrames = vnStimulusSeqID == nThisStimID;
   vbThisStimFrames = vbThisStimFrames & ...
      (vtTimeInStimPresentation >= mtUseStimDurations(nStim, 1)) & ...
      (vtTimeInStimPresentation <= mtUseStimDurations(nStim, 2));
   
   % - Extract response and compute mean
   mfResponseMean = double(oStack.SummedAlignedFrames(:, :, vbThisStimFrames, nChannel)) ./ nnz(vbThisStimFrames);

   % - Get blank frames for this stimuli
   [tfBlankMeans, tfBlankStds, vnBlankFrameIndices] = ...
      GetCorrespondingBlankFrames(oStack, {':' ':' vbThisStimFrames});
   
   % - Compute corrected std.dev. for combined blank
   tfBlankMeanSqr = tfBlankMeans.^2;
   tfBlankCorrPart = tfBlankMeanSqr + tfBlankStds.^2;
   clear vfWeights;
   vfWeights(1, 1, :) = accumarray(vnBlankFrameIndices, 1);
   vnPixelSize = size(tfBlankMeanSqr);
   mfCorrStdBlank = real(sqrt(sum(tfBlankCorrPart .* repmat(vfWeights, [vnPixelSize([1 2]) 1]), 3) / sum(vfWeights) - mean(tfBlankMeanSqr, 3)));
   
   % - Compute mean blank-subtracted response
   mfStimMeanBlankSub = mfResponseMean - mean(tfBlankMeans, 3);
   
   % - Z score is (mean obs. - mean blank) / (std blank / sqrt(num stim frames))
   %   This uses the expect standard error of the mean of several response
   %   frames
   tfStimZ(:, :, nStim) = mfStimMeanBlankSub ./ (mfCorrStdBlank / sqrt(nnz(vbThisStimFrames))); %#ok<AGROW,SNASGU>
end

% -- Spatially filter Z score matrices

fBPFiltWidth = 5; % pixels

hLowPass = fspecial('gaussian', fBPFiltWidth*3, fBPFiltWidth/4);

tfStimZ(isnan(tfStimZ) | isinf(tfStimZ)) = 0;
tfStimZFilt = imfilter(tfStimZ, hLowPass);

% - Find maximum Z score over all stimulus responses
mfMaxStimZ = max(tfStimZ .* (tfStimZ>=0), [], 3);
% mfMinStimZ = min(tfStimZ .* (tfStimZ<=0), [], 3);

mfMaxStimZFilt = max(tfStimZFilt .* (tfStimZ>=0), [], 3);
% mfMinStimZFilt = min(tfStimZFilt .* (tfStimZ<=0), [], 3);


% -- Identify responsive pixels
% - Responsive pixels are those for whom the stimulus Z-score passes an alpha threshold

fMaxZScoreThreshold = norminv(1-fAlphaThreshold, 0, 1);
% fMinZScoreThreshold = -fMaxZScoreThreshold;              % Symmetry

% - Ignore pixels removed by the alignment mask
mbAlignMask = oStack.GetAlignedMask;

mbRespMask = (mfMaxStimZFilt > fMaxZScoreThreshold) & mbAlignMask;% | (mfMinStimZFilt < fMinZScoreThreshold);

sRespRegions = bwconncomp(mbRespMask);


%% -- Group pixels based on response pattern

% mbStimResponse = reshape(tfStimZFilt > fMaxZScoreThreshold, prod(vnStackSize(1:2)), []);
% 
% % - Create a new region structure
% sFiltRespRegions.Connectivity = 8;
% sFiltRespRegions.ImageSize = sRespRegions.ImageSize;
% sFiltRespRegions.NumObjects = 0;
% sFiltRespRegions.PixelIdxList = {};
% 
% % - Examine each region in turn
% for (nRegion = sRespRegions.NumObjects:-1:1)
%    % - Extract response patterns for this region
%    mbResponsePatterns = mbStimResponse(sRespRegions.PixelIdxList{nRegion}, :);
%    
%    % - Find unique response patterns
%    [mbUniquePatterns, nul, vnMatchingIndices] = unique(mbResponsePatterns, 'rows');
%    
%    % - Loop over patterns and assign to new objects
%    for (nPattern = 1:size(mbUniquePatterns, 1))
%       % - Find pixels that respond with this pattern
%       vbMatchingPattern = vnMatchingIndices == nPattern;
%       vnMatchingPixels = sRespRegions.PixelIdxList{nRegion}(vbMatchingPattern);
%       
%       % - Make a mask including only these pixels
%       mbTestMask = false(sFiltRespRegions.ImageSize);
%       mbTestMask(vnMatchingPixels) = true;
%       
%       % - Split these pixels into regions and record
%       sTheseRegions = bwconncomp(mbTestMask);
%       sFiltRespRegions.PixelIdxList = [sFiltRespRegions.PixelIdxList(:)' sTheseRegions.PixelIdxList(:)'];
%    end
% end
% 
% sRespRegions = sFiltRespRegions;
% sRespRegions.NumObjects = numel(sRespRegions.PixelIdxList);


%% -- Filter regions by size


if (exist('nMinRegionSize', 'var') && ~isempty(nMinRegionSize))
   vnRegionSizes = cellfun(@numel, sRespRegions.PixelIdxList);
   vbAcceptRegion = vnRegionSizes >= nMinRegionSize;
   sRespRegions.PixelIdxList = sRespRegions.PixelIdxList(vbAcceptRegion);
   sRespRegions.NumObjects = numel(sRespRegions.PixelIdxList);
end

if (exist('nMaxRegionSize', 'var') && ~isempty(nMaxRegionSize))
   vnRegionSizes = cellfun(@numel, sRespRegions.PixelIdxList);
   vbAcceptRegion = vnRegionSizes <= nMaxRegionSize;
   sRespRegions.PixelIdxList = sRespRegions.PixelIdxList(vbAcceptRegion);
   sRespRegions.NumObjects = numel(sRespRegions.PixelIdxList);
end

% - Restore mask matrix
mbRespMask = labelmatrix(sRespRegions) > 0;


% -- Order regions by Z-score

% - Find max Z score for each region
vfRegionZScore = cellfun(@(r)(max(mfMaxStimZFilt(r))), sRespRegions.PixelIdxList);

% - Convert to alphas
vfRegionAlpha = 1-normcdf(vfRegionZScore);

% - Sort regions by max Z score
[nul, vnSortOrder] = sort(abs(vfRegionZScore), 'descend');
sRespRegions.PixelIdxList = sRespRegions.PixelIdxList(vnSortOrder);
sRespRegions.vfRegionZScore = vfRegionZScore(vnSortOrder);
sRespRegions.vfRegionAlpha = vfRegionAlpha(vnSortOrder);

% - Restore DFF conversion to stack
oStack.BlankNormalisation(strNormalisation);


% --- END of ResponsiveMask.m ---
