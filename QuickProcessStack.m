function QuickProcessStack(strFilename)

% QuickProcessStack - FUNCTION Quickly average and process a set of focus stacks
%
% Usage: QuickProcessStack(strFilename)
%        QuickProcessStack('config')
%
% This function is designed to act as a drag-and-drop target along with the
% helper functions 'openfcs', 'opentif' and 'opentiff'.  It will
% automatically process a set of dropped FocusStack files, de-randomise
% stimulus order (if required) and display trial-averaged responses.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 28th October, 2011

% -- Shared variables

global QPS_sConfig;
persistent cellFilenames tTimer;


% -- Check arguments

if (nargin == 0)
   disp('*** QuickProcessStack: Incorrect usage');
   help QuickProcessStack;
   return;
end

strFilename = CellFlatten({strFilename});

switch lower(strFilename{1})
   case {'config', 'configure', 'configuration'}
      QuickProcessStack_config;
      
   otherwise
      % - Append this filename(s)
      cellFilenames = [cellFilenames strFilename(:)];
      
      % - Either create or extend the timer
      if isempty(tTimer)
         tTimer = timer('Name', 'FocusStack open execution timer', 'StartDelay', 0.5, 'TimerFcn', @(o,e)QuickProcessStack_openfiles);
         start(tTimer);
      else
         stop(tTimer);
         start(tTimer);
      end
end


% --- END of QuickProcessStack FUNCTION ---

   function QuickProcessStack_openfiles
      % - Does the configuration exist?
      while (isempty(QPS_sConfig))
         uiwait(QuickProcessStack_config);
      end      
      
      % - Remove timer
      stop(tTimer);
      delete(tTimer);
      tTimer = [];
      
      % - Try to create a focus stack
      try
         disp('--- QuickProcessStack: Creating FocusStack...');
         fsStack = FocusStack(cellFilenames);
         
         % - Clear filenames
         cellFilenames = [];
         
         % - Assign stack annotations, if required
         if (isempty(fsStack.tFrameDuration))
            fsStack.tFrameDuration = QPS_sConfig.tFrameDuration;
         end
         
         if (all(cellfun(@(c)all(isnan(c)), fsStack.cvnSequenceIDs)))
            vnSequenceIDs = [1:QPS_sConfig.nNumStimuli nan]; 
            fsStack.cvnSequenceIDs = repmat({vnSequenceIDs}, 1, numel(fsStack.cstrFilenames));
         end
         
         if (isempty(fsStack.vtStimulusDurations))
            fsStack.vtStimulusDurations = QPS_sConfig.vtStimulusDurations;
         end
         
         % - Always assign use times to override the stack default (which is to take all frames)
         fsStack.mtStimulusUseTimes = QPS_sConfig.mtStimulusUseTimes;
         
         % - Align stack
         if (QPS_sConfig.bAlign)
            disp('--- QuickProcessStack: Aligning stack...');
            fsStack.Align(1, false, 10);
         end
         
         % - Assign black trace
         if (QPS_sConfig.bAssignBlack)
            disp('--- QuickProcessStack: Assigning black trace...');
            fsStack.DefineBlackRegion;
         end
         
         % - Assign blank for blank subtraction
         if (QPS_sConfig.bAssignBlank)
            disp('--- QuickProcessStack: Assigning blank frames...');
            AssignBlank(fsStack);
         end
         
         % - Are there regions defined?
         if (isfield(QPS_sConfig, 'sRegions') && ~isempty(QPS_sConfig.sRegions))
            sRegions = QPS_sConfig.sRegions;
         else
            % - Automatically define regions
            disp('--- QuickProcessStack: Identifying ROIs...');
            if (size(fsStack, 4) > 1)
               sRegions = FindCells_GRChannels(fsStack);
            else
               sRegions = FindCells_GChannel(fsStack);
            end
         end
         
         % - Add a "neuropil" region
         nNumPixels = prod(size(fsStack, 1:2)); %#ok<PSIZE>
         sRegionsPlusNP = sRegions;
         sRegionsPlusNP.NumObjects = sRegionsPlusNP.NumObjects+1;
         sRegionsPlusNP.PixelIdxList = [{1:nNumPixels} sRegionsPlusNP.PixelIdxList];
         
         % - Extract responses, order by Z-score
         disp('--- QuickProcessStack: Extracting responses...');
         
         if (QPS_sConfig.bAssignBlank)
            fsStack.BlankNormalisation('subtractive');
         end
         
         [vfBlankStds, mfStimMeanResponses, mfStimStds, mfRegionTraces, tfTrialResponses, tnTrialSampleSizes] = ...
            ExtractRegionResponses(fsStack, sRegionsPlusNP, QPS_sConfig.nBlankStimID, QPS_sConfig.fhExtract);
         
         if (QPS_sConfig.bAssignBlank)
            % - Estimate neuropil-rejecting Z score thresholds
            disp('--- QuickProcessStack: Estimating neuropil distribution...');
            [sNeuropilZThresh] = ...
               ThresholdNeuropilResponse(fsStack, sRegions, QPS_sConfig.vnUseStimIDs, 0.01, 1, QPS_sConfig.fhExtract);
            
            % - Compute Z-score against blank STDs
            mnStimSampleSizes = sum(tnTrialSampleSizes, 3);
            % vnROINumPixels = cellfun(@numel, sRegionsPlusNP.PixelIdxList)';
            mfBlankStdsCorr = repmat(vfBlankStds', 1, QPS_sConfig.nNumStimuli) ./ sqrt(mnStimSampleSizes);
            % mfStimZScores = mfStimMeanResponses ./ mfBlankStdsCorr ./ repmat(sqrt(vnROINumPixels), 1, QPS_sConfig.nNumStimuli);
            mfStimZScores = mfStimMeanResponses ./ mfBlankStdsCorr;
            
            tfBlankStdsCorr = repmat(vfBlankStds', [1, QPS_sConfig.nNumStimuli, size(tfTrialResponses, 3)]) ./ sqrt(tnTrialSampleSizes);
            % tfStimZScoresTrials = tfTrialResponses ./ tfBlankStdsCorr ./ repmat(sqrt(vnROINumPixels), 1, QPS_sConfig.nNumStimuli);
            tfStimZScoresTrials = tfTrialResponses ./ tfBlankStdsCorr;
            
            % - Measure significance of response (above neuropil)
            nNumTrials = numel(fsStack.cstrFilenames);
            mbReliableResponse = sum(tfStimZScoresTrials > sNeuropilZThresh.fTrialStimResp, 3) > (nNumTrials/2);
            vbSignificantResponse = nanmax(mfStimZScores, [], 2) > sNeuropilZThresh.fMaxStimResp;
            vbKeepROI = vbSignificantResponse & any(mbReliableResponse, 2);
            
            % - Include NP
            vbKeepROI(1) = true;
            
            % - Sort by Z-score, removing unresposive cells
            vfMaxZScores = max(mfStimZScores, [], 2);
            vfMaxZScores(1) = inf;  % - Force NP inclusion, as first region
            [nul, vnSort] = sort(vfMaxZScores.*vbKeepROI, 'descend');
            
            sRegionsPlusNP.PixelIdxList = sRegionsPlusNP.PixelIdxList(vnSort);            
            mfRegionTraces = mfRegionTraces(vnSort, :);
            tfTrialResponses = tfTrialResponses(vnSort, :, :); %#ok<NASGU>
            mbReliableResponse = mbReliableResponse(vnSort, :);
            vbKeepROI = vbKeepROI(vnSort);
            
         else
            vbKeepROI = true(1, sRegionsPlusNP.NumObjects);
            mbReliableResponse = [];
         end
         
         disp('--- QuickProcessStack: Plotting figures...');
         
         % - Make an overview plot
         vnPlotROIs = find(vbKeepROI);
         vhFigures = ...
            MakeROIOverviewFigure(fsStack, [], 1:size(fsStack, 4), sRegionsPlusNP, vnPlotROIs(2:end), true);
         
         % - Plot region responses
         vhFigures = [vhFigures ...
            PlotRegionResponses(mfRegionTraces, fsStack, sRegionsPlusNP, ...
                             QPS_sConfig.mtStimulusLabelTimes, 0.1, 5, ...
                             nnz(vbKeepROI), [], mbReliableResponse)];

         % - Assign stack to base workspace
         assignin('base', 'fsQuickStack', fsStack);
         
         % - Tile figures
%          TileFigures(vhFigures);
         
         disp('--- QuickProcessStack: Done.');
                          
      catch mErr
         % - Clean up
         cellFilenames = [];
         rethrow(mErr);
      end      
   end
end

function AssignBlank(fsStack)

global QPS_sConfig;

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
   % - Get the average blank frame for this block
   vbBlockFrames = vnBlockIndex == nBlock;
   vbBlockBlankStimFrames = vbBlockFrames & (vnStimulusSeqID == QPS_sConfig.nBlankStimID);
   vbBlockBlankFrames = vbBlockBlankStimFrames & vbUseFrame;
   
   tfBlankFrames = double(fsStack.AlignedStack(:, :, vbBlockBlankFrames, 1));
   mfThisBlankAvg = nanmean(tfBlankFrames, 3);
   
   % - Filter zeros, replace with NaN
   mfThisBlankAvg(mfThisBlankAvg == 0) = nan;
   
   % - Calculate std. dev. for divisive normalisation
   mfThisBlankStd = nanstd(tfBlankFrames ./ repmat(mfThisBlankAvg, [1 1 nnz(vbBlockBlankFrames)]), [], 3);
   
   % - Assign the blank frame to the whole block (by default)
   fsStack.AssignBlankFrame(cat(3, mfThisBlankAvg, mfThisBlankStd), vbBlockFrames);
   
   % - Assign separate blank frames for each other stimulus segment
   for (nStimSeqID = QPS_sConfig.vnUseStimIDs(:)')
      % - Find all frames for this presentation
      vbPresFrames = vbBlockFrames & (vnStimulusSeqID == nStimSeqID);
      
      % - Find blank frames for this presentation
      vbPresBlankFrames =  vbPresFrames & ...
         (vtTimeInStimPresentation >= QPS_sConfig.mtBlankTimes(nStimSeqID, 1)) & ...
         (vtTimeInStimPresentation <= QPS_sConfig.mtBlankTimes(nStimSeqID, 2));
      
      % - Extract these blank frames
      tfPresBlank = fsStack.AlignedStack(:, :, vbPresBlankFrames, 1);
      
      % - Assign only the mean; Std.Dev is taken from block blank
      fsStack.AssignBlankFrame(cat(3, nanmean(tfPresBlank, 3), mfThisBlankStd), vbPresFrames);
   end
end

% - Restore warnings
warning(wOld);

end

function TileFigures(vhFigures)
    
    N = numel(vhFigures);
    
    C = ceil( sqrt(N) );
    R = ceil( N/C );
    set( vhFigures, 'units', 'normalized' );
    % set( hFig, 'toolbar', 'none' );    
    % set( hFig, 'toolbar', 'figure' );
    w = 1/C;
    h = 0.99/R;
    
    for n = 1 : N
        [r c] = ind2sub( [R C], n );
        p = [ (c-1)*w 0.01+(r-1)*h,  w*1, h*0.82 ];
        set( vhFigures(n),'position',p );     
        figure( vhFigures(n) );
        
    end
end

% --- END QuickProcessStack.m ---
