function [vfBlankStds, mfStimMeanResponses, mfStimStds, ...
          mfRegionTraces, tfTrialResponses, tnFramesInSample, cvfTrialTraces] = ...
   ExtractRegionResponses(fsStack, sRegions, nBlankStimID, fhExtractionFunction, tBaseTimeShift)

% ExtractRegionResponses - FUNCTION Extract responses from identified regions of interest
%
% Usage: [vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, ...
%         tfTrialResponses, tnFramesInSample, cvfTrialTraces] = ...
%           ExtractRegionResponses(fsStack, sRegions, nBlankStimID, <, fhExtractionFunction, tBaseTimeShift>)
%
% 'fsStack' is a FocusStack object.  If stimulus information is embedded (ie
% 'vtStimulusDirations', 'vtStimulusStartTimes', 'vtStimulusEndTimes',
% 'mtStimulusUseTimes') then the time course of the stack will be segmented by
% stimuli.  This segmentation applies to the return variables 'mfStimMeans',
% 'mfStimStds', 'mnStimSampleSizes', 'tfTrialResponses'.
%
% 'sRegions' is a ROI structure, as returned from bwconncomp.  Information will
% be returned for each region.
%
% The optional argument 'fhExtractFunction' is a function handle that extract
% the desired components of the response trace, and condenses the response to a
% single number.  It also determines how the channels in 'fsStack' are handled.
% The function must have the signature:
% [mfRawTrace, vfRegionTrace, fResponse, nFramesInSample] = @fh(fsData, vnPixels, vnFrames)
% The function 'fh' must accept a FocusStack object, the linear index of the
% pixels that define an ROI, and the frame indices that define which frames are
% included in a response.  The function must reference 'fsStack', also
% determining which channels of the stack must be used, and return a scalar
% description of the response from those frames.  For example, a function might
% compute the average over time and space of the response trace of channel 1.
% Or the function might compute the ratio of channels 1 and 2, find the peak of
% the response, and average a few points around the peak.
% 'vfTrace' must be the response trace over time for the defined region, as
% considered by the function.  In the default case, this is the space-averaged
% data from channel 1.
%
% If 'fhExtractionFunction' is not supplied, the default will be to extract the
% time and space average of channel 1.
%
% The optional parameter 'tBaseTimeShift' shifts the effective stimulus
% start and end times by the provided amount (in seconds; positive or
% negative).  Use this to shif the analysed frames in the stack earlier or
% later. (Default: Not used).  See FrameStimulusInfo for more information.
%
% 'vfBlankStds' will be a vector [Rx1], where 'R' is the number of ROIs.  Each
% element contains the blank standard deviation (corrected for the number of
% pixels in each ROI) PER FRAME for a given ROI.  This value must be further
% corrected to compare with averages over several frames (ie divide by
% sqrt(nNumFrames)).
%
% 'mfStimMeanResponses' will be a matrix [RxS], where 'R' is the number of ROIs
% and 'S' is the number of stimuli.  Each element in 'mfStimResponses' is the
% mean of the scalar values returned by the extraction function
% 'fhExtractionFunction', over all trials.  In the default case, this will be
% the space-time mean of all pixels within each ROI, during the stimulus
% presentation time, extracted from channel 1 only of the stack.
%
% 'mfStimStds' will be a matrix [RxS], where each element corresponds to
% 'mfStimMeanResponses'.  It will contain the standard deviation of the trial
% responses for each ROI and each stimulus.
%
% 'mfRegionTraces' will be a matrix [RxT], where each row of 'R' (the number of
% ROIs) contains the response trace over time for a single region.  'T' is the
% total number of frames in the stack.  These responses will be as extracted
% by 'fhExtractionFunction'.
%
% 'tfTrialResponses' will be a matrix [RxSxM], where 'R' is the number of ROIs,
% 'S' is the number of stimuli and 'M' is the number of trials.  Each element of
% 'tfTrialResponses' is the scalar returned by 'fhExtractionFunction' for a
% single ROI, for a single stimulus, for a single trial.
%
% 'tnFramesInSample' will be a matrix [RxSxM], conventions as above, where
% each element contains the number of frames used to form the corresponding
% element in 'tfTrialResponses'.  For example, it could be the number of
% samples averaged together to form a mean.
%
% 'cvfTrialTraces' will be a cell tensor {RxSxM}, conventions as above.
% Each cell element will contain the raw data trace for the corresponding
% region, stimulus and trial.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 23rd november, 2010

% -- Check arguments

if (nargin < 2)
   disp('*** ExtractRegionResponses: Incorrect usage');
   help ExtractRegionResponses;
end

if (~exist('fhExtractionFunction', 'var') || isempty(fhExtractionFunction))
   fhExtractionFunction = ExtractMean;
end

if (~exist('tBaseTimeShift', 'var'))
   tBaseTimeShift = [];
end

nNumFrames = size(fsStack, 3);

% - Check that required stack data are present
if (  isempty(fsStack.cvnSequenceIDs) || ...
      isempty(fsStack.vtStimulusDurations) || ...
      isempty(fsStack.tFrameDuration) || ... 
      isempty(fsStack.mtStimulusUseTimes))
   disp('--- ExtractRegionResponses: Warning: Not all required stack data is available.');
   bSegmentStack = false;
else
   
   bSegmentStack = true;
   nNumStimuli = fsStack.nNumStimuli;
   [vtGlobalTime, ...
      vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
      vnStimulusSeqID, vtTimeInStimPresentation, ...
      vnPresentationIndex, vbUseFrame] = ...
         FrameStimulusInfo(fsStack, 1:nNumFrames, tBaseTimeShift);
end


% -- Average region traces together and extract regions

nNumRegions = sRegions.NumObjects;
fprintf('ExtractRegionResponses...');

% - Extract the region traces
nNumFrames = size(fsStack, 3);
[nul, cvfRegionTrace] = fhExtractionFunction(fsStack, sRegions.PixelIdxList, 1:nNumFrames);
if ~iscell(cvfRegionTrace)
    mfRegionTraces = cvfRegionTrace;
else
    mfRegionTraces = vertcat(cvfRegionTrace{:});
end

% - Extract the region blanks
if (~isempty(nBlankStimID))
   vbBlankFrames = (vnStimulusSeqID == nBlankStimID) & vbUseFrame;
else
   vbBlankFrames = ~vbUseFrame;
end

[cmfBlankTrace, cvfRegionBlank] = fhExtractionFunction(fsStack, sRegions.PixelIdxList, vbBlankFrames);
if ~iscell(cmfBlankTrace)
    cmfBlankTrace = {cmfBlankTrace};
    cvfRegionBlank = {cvfRegionBlank};
end

vfBlankStds = cellfun(@nanstd, cvfRegionBlank);


% -- Extract region responses for each stimulus and trial

if (bSegmentStack)
   % -- Extract segments for each stimulus and average
   fprintf('\b\b\b: %4d%%...', 0);

   % - Allocate single-trial storage
   nNumBlocks = numel(fsStack.cstrFilenames);
   tfTrialResponses = nan(nNumRegions, nNumStimuli, nNumBlocks);
   cvfTrialTraces = cell(nNumRegions, nNumStimuli, nNumBlocks);
   
   for (nStimSeqID = nNumStimuli:-1:1)
      % - Peel off trials
      for (nTrialID = nNumBlocks:-1:1)
         % - Find matching stack frames for this trial
         vbFramesThisTrialStim = (vnStimulusSeqID == nStimSeqID) & (nTrialID == vnBlockIndex) & vbUseFrame;
         [nul, cvfTrialTracesAll, cfTrialResponses, cnFramesInSample] = ...
            fhExtractionFunction(fsStack, sRegions.PixelIdxList, vbFramesThisTrialStim);

        if ~iscell(cvfTrialTracesAll)
            cvfTrialTracesAll = {cvfTrialTracesAll};
            cfTrialResponses = {cfTrialResponses};
            cnFramesInSample = {cnFramesInSample};
        end
        
         % - Separate out regions
         [cvfTrialTraces{:, nStimSeqID, nTrialID}] = deal(cvfTrialTracesAll{:});
         tfTrialResponses(:, nStimSeqID, nTrialID) = [cfTrialResponses{:}];
         tnFramesInSample(:, nStimSeqID, nTrialID) = [cnFramesInSample{:}];
      end
      fprintf('\b\b\b\b\b\b\b\b%4d%%...', round((nNumStimuli - nStimSeqID+1) / nNumStimuli * 100));
   end
   fprintf('\b\b\b\b\b\b\b\b%4d%%.\n', 100);
   drawnow;
end


if (bSegmentStack)
   % - Calculate averages and std devs
   mfStimMeanResponses = mean(tfTrialResponses, 3);
   mfStimStds = std(tfTrialResponses, [], 3);

else
   mfStimMeanResponses = [];
   mfStimStds = [];
   tfTrialResponses = [];
end

% --- END of ExtractRegionResponse.m ---
