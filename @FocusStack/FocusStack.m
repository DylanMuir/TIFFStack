% FocusStack - CLASSDEF Contains memory-mapped two-photon stacks on disk
%
% To construct a FocusStack object:
%    fsStack = FocusStack(cstrFilenames, bWritable)

classdef FocusStack < handle
   properties (SetAccess = private)
      % bWritable - Can we write back to the file no disk?
      bWritable = false;
      cstrFilenames = {};
      bConvertToDFF = false;
      bSubtractBlank = false;
   end
   
   properties (SetAccess = private, GetAccess = private)
      vnNumFrames;
      nNumChannels;
      vnFrameSize = [];
      strDataClass;
   end
   
   properties (SetAccess = private, GetAccess = private, Transient = true)
      vhMemMapFileHandles;
      vbCachedAlignedFrames = [];
      oAlignedFrameCache = [];
      cmfBlankFrames = {};
      mnAssignedBlankMeanFrames = [];
      mnAssignedBlankStdFrames = [];
      bAssignedSequenceIDs = false;
      bAssignedStimStartTimes = false;
      bAssignedStimEndTimes = false;
      vsHeaders;
   end

   properties (Dependent = true, SetAccess = private)
      vnStimulusIDs;
      tBlankTime;
      nNumStimuli;
   end

   properties
      fPixelsPerUM = [];
      tFrameDuration = [];
      fZStep = [];
      cvnSequenceIDs = {};
      vtStimulusDurations = [];
      vtStimulusStartTimes = [];
      vtStimulusEndTimes = [];
      mtStimulusUseTimes = [];
      mfFrameShifts;
      vfBlackTrace;
      bSubtractBlack = false;
   end
   
%% --- External private method declarations
   
   methods (Access = private)
      OpenFiles(oStack)
      CloseFiles(oStack)
      [cvnFileRefs, cvnFullRefs, vnDataSize] = GetFullFileRefs(oStack, cRefs)
   end

%% --- Constructor and destructor
   
   methods
      % FocusStack - CONSTRUCTOR
      function oStack = FocusStack(cstrFilenames, bWritable)
         % - Assign default values if no arguments were provided
         if (nargin == 0)
            return;
         end
         
         % - Assign file names
         if (~iscellstr(cstrFilenames))
            oStack.cstrFilenames = {cstrFilenames};
         else
            oStack.cstrFilenames = cstrFilenames;
         end

         % - Assign 'writable' flag
         if (~exist('bWritable', 'var') || isempty(bWritable))
            oStack.bWritable = false;
         end
         
         % - Try to open the files
         OpenFiles(oStack);
      end
      
      % delete - DESTRUCTOR
      function delete(oStack)
         % - Close the file handles
         CloseFiles(oStack);
      end
      
%% --- Property accessors
      
      % - BlankNormalisation - METHOD Set/Get blank normalisation mode
      %
      % Usage: [strOldNormalisation] = BlankNormalisation(oStack, strNormalisation)
      function [strOldNormalisation] = BlankNormalisation(oStack, strNormalisation)
         if (oStack.bConvertToDFF)
            strOldNormalisation = 'divisive';
         
         elseif (oStack.bSubtractBlank)
            strOldNormalisation = 'subtractive';
            
         else
            strOldNormalisation = 'none';
         end
         
         if (exist('strNormalisation', 'var'))
            switch(strNormalisation)
               case {'divisive', 'divide'}
                  oStack.bConvertToDFF = true;
                  oStacl.bSubtractBlank = false;
                  
               case {'subtractive', 'subtract'}
                  oStack.bConvertToDFF = false;
                  oStack.bSubtractBlank = true;
                  
               case 'none'
                  oStack.bConvertToDFF = false;
                  oStack.bSubtractBlank = false;
                  
               otherwise
                  error('FocusStack:BadArgument', ...
                     '*** FocusStack/BlankNormalisation: Normalisation type must be one of {''divisive'', ''subtractive'', ''none''}.');
            end
         end
         
         % - Clear return argument, if not requested
         if (nargout == 0)
            clear strOldNormalisation;
         end
      end

      % set.mfFrameShifts - SETTER for 'mfFrameShifts'
      function set.mfFrameShifts(oStack, mfFrameShifts)
         % -- Check data
         if (~isnumeric(mfFrameShifts))
            error('FocusStack:InvalidFrameShiftData', '*** FocusStack/set.mfFrameShifts: The assigned frame shifts must be numeric.');
         end
         
         if (~isequal(size(mfFrameShifts), [sum(oStack.vnNumFrames) 2])) %#ok<MCSUP>
            error('FocusStack:InvalidFrameShiftData', ...
               '*** FocusStack/set.mfFrameShifts: The assigned frame shifts must be size [%d 2] for this stack.', sum(oStack.vnNumFrames)); %#ok<MCSUP>
         end
         
         % - Assign data
         oStack.mfFrameShifts = mfFrameShifts;
         
         % - Reset alignment cache
         oStack.oAlignedFrameCache = []; %#ok<MCSUP>
      end
      
      % set.vfBlackTrace - SETTER for 'vfBlackTrace'
      function set.vfBlackTrace(oStack, vfBlackTrace)
         % -- Check data
         if (~isnumeric(vfBlackTrace))
            error('FocusStack:InvalidBlackTraceData', ...
               '*** FocusStack/set.vfBlackTrace: The assigned black trace must be numeric.');
         end
         
         if (numel(vfBlackTrace) ~= size(oStack, 3))
            error('FocusStack:InvalidBlackTraceData', ...
               '*** FocusStack/set.vfBlackTrace: The assigned black trace must have [%d] elements for this stack.', size(oStack, 3));
         end
         
         % - Assign data
         oStack.vfBlackTrace = double(reshape(vfBlackTrace, [], 1));
         
         % - Switch on black subtraction
         oStack.bSubtractBlack = true; %#ok<MCSUP>
      end
      
      % get.vnStimulusIDs - GETTER for 'vnStimulusIDs'
      function [vnStimulusIDs] = get.vnStimulusIDs(oStack)
         vnStimulusIDs = [oStack.vsHeaders.nStimulusID];
      end
      
      % get.fBlankTime - GETTER for 'tBlankTime'
      function [tBlankTime] = get.tBlankTime(oStack)
         vtBlankTimes = [oStack.vsHeaders.tBlankTime];
         
         if (all(vtBlankTimes == vtBlankTimes(1)))
            tBlankTime = vtBlankTimes(1);
         else
            tBlankTime = vtBlankTimes;
         end
      end
      
      % set.cvnSequenceIDs - SETTER for 'cvnSequenceIDs'
      function set.cvnSequenceIDs(oStack, cvnSequenceIDs)
         if (  ~iscell(cvnSequenceIDs) || ...
               numel(cvnSequenceIDs) ~= numel(oStack.cstrFilenames)) %#ok<MCSUP>
            error('FocusStack:InvalidArgument', ...
               '*** FocusStack/set.cvnSequenceIDs: ''cvnSequenceIDs'' must be a cell array with [%d] element(s) for this stack.', ...
               numel(oStack.cstrFilenames)); %#ok<MCSUP>
         end
         
         % - Reshape all sequence ID arrays
         cvnSequenceIDs = cellfun(@(c)(reshape(c, [], 1)), cvnSequenceIDs, 'UniformOutput', false);
         
         % - Assign cvnSequenceIDs
         oStack.cvnSequenceIDs = reshape(cvnSequenceIDs, [], 1);
         oStack.bAssignedSequenceIDs = true; %#ok<MCSUP>
      end        
      
      % get.cvnSequenceIDs - GETTER for 'cvnSequenceIDs'
      function [cvnSequenceIDs] = get.cvnSequenceIDs(oStack)
         % - Have the sequence IDs been assigned?
         if (oStack.bAssignedSequenceIDs)
            cvnSequenceIDs = oStack.cvnSequenceIDs;
            
         else
            % - Compute them from the block headers
            cvnSequenceIDs = {oStack.vsHeaders.vnSequenceIDs};
            
            % - Add a "skip to end of block" tag to the end of each block
            cvnSequenceIDs = cellfun(@(c)([reshape(c, [], 1); nan]), cvnSequenceIDs, 'UniformOutput', false);
         end
      end
      
      % set.bSubtractBlack - SETTER for 'bSubtractBlack'
      function set.bSubtractBlack(oStack, bSubtractBlack)
         if (~islogical(bSubtractBlack) || ~isscalar(bSubtractBlack))
            error('FocusStack:InvalidArgument', ...
               '***FocusStack/set.bSubtractBlack: ''bSubtractBlack'' must be a scalar logical.');
         end
         
         oStack.bSubtractBlack = bSubtractBlack;
      end
      
      % set.bConvertToDFF - SETTER for 'bConvertToDFF'
      function set.bConvertToDFF(oStack, bConvertToDFF)
         if (~islogical(bConvertToDFF) || ~isscalar(bConvertToDFF))
            error('FocusStack:InvalidArgument', ... 
               '***FocusStack/set.bConvertToDFF: ''bConvertToDFF'' must be a scalar logical.');
         end
         
         % Assign bConvertToDFF
         oStack.bConvertToDFF = bConvertToDFF;
      end
      
      % get.nNumStimuli - GETTER for 'nNumStimuli'
      function [nNumStimuli] = get.nNumStimuli(oStack)
         vnSequenceIDs = [oStack.cvnSequenceIDs{:}];
         nNumStimuli = numel(unique(vnSequenceIDs(~isnan(vnSequenceIDs))));
      end
      
      % set.vtStimulusDurations - SETTER for 'vtStimulusDurations'
      function set.vtStimulusDurations(oStack, vtStimulusDurations)
         if (~isnumeric(vtStimulusDurations) || (numel(vtStimulusDurations) ~= oStack.nNumStimuli)) %#ok<MCSUP>
            error('FocusStack:InvalidArgument', ...
               '*** FocusStack/set.vtStimulusDurations: ''vtStimulusDurations'' must have [%d] elements for this stack.', ...
               oStack.nNumStimuli); %#ok<MCSUP>
         end
         
         % - Assign vtStimulusDurations
         oStack.vtStimulusDurations = reshape(vtStimulusDurations, [], 1);
      end
      
      % set.vtStimulusStartTimes - SETTER for 'vtStimulusStartTimes'
      function set.vtStimulusStartTimes(oStack, vtStimulusStartTimes)
         cvnSequenceIDs = oStack.cvnSequenceIDs; %#ok<MCSUP,PROP>
         vnStimOrder = vertcat(cvnSequenceIDs{:}); %#ok<PROP>
         nNumPresentations = numel(vnStimOrder);
         
         if (~isnumeric(vtStimulusStartTimes) || (numel(vtStimulusStartTimes) ~= nNumPresentations) && ~isempty(vtStimulusStartTimes)) 
            error('FocusStack:InvalidArgument', ...
               '*** FocusStack/set.vtStimulusStartTimes: ''vtStimulusStartTimes'' must have [%d] elements for this stack.', ...
               oStack.nNumStimuli); %#ok<MCSUP>
         end
         
         % - Assign vtStimulusStartTimes
         oStack.vtStimulusStartTimes = reshape(vtStimulusStartTimes, [], 1);
         oStack.bAssignedStimStartTimes = ~isempty(vtStimulusStartTimes); %#ok<MCSUP>
      end
      
      % get.vtStimulusStartTimes - GETTER for 'vtStimulusStartTimes'
      function [vtStimulusStartTimes] = get.vtStimulusStartTimes(oStack)

         % - If the stimulus start times have been assigned, just return them
         if (oStack.bAssignedStimStartTimes)
            vtStimulusStartTimes = oStack.vtStimulusStartTimes;
         
         else
            % - If the stimulus start times are not already computed, assume that
            % each stimulus follows the other with no gap.
            vtStimulusStartTimes = calc_stim_start_end_times(oStack);
         end
      end
      
      % set.vtStimulusEndTimes - SETTER for 'vtStimulusEndTimes'
      function set.vtStimulusEndTimes(oStack, vtStimulusEndTimes)
         cvnSequenceIDs = oStack.cvnSequenceIDs; %#ok<MCSUP,PROP>
         vnStimOrder = vertcat(cvnSequenceIDs{:}); %#ok<PROP>
         nNumPresentations = numel(vnStimOrder);

         if (~isnumeric(vtStimulusEndTimes) || (numel(vtStimulusEndTimes) ~= nNumPresentations) && ~isempty(vtStimulusEndTimes)) %#ok<MCSUP>
            error('FocusStack:InvalidArgument', ...
               '*** FocusStack/set.vtStimulusEndTimes: ''vtStimulusEndTimes'' must have [%d] elements for this stack.', ...
               oStack.nNumStimuli); %#ok<MCSUP>
         end
         
         % - Check values for end times
         if (~isempty(oStack.vtStimulusStartTimes) && ~isempty(vtStimulusEndTimes)) %#ok<MCSUP>
            if (any(vtStimulusEndTimes < oStack.vtStimulusStartTimes)) %#ok<MCSUP>
               error('FocusStack:InvalidArgument', ...
                  '*** FocusStack/set.vtStimulusEndTimes: ''vtStimulusEndTimes'' must all fall after the corresponding start time.');
            end
         end
         
         % - Assign vtStimulusStartTimes
         oStack.vtStimulusEndTimes = reshape(vtStimulusEndTimes, [], 1);
         oStack.bAssignedStimEndTimes = ~isempty(vtStimulusEndTimes); %#ok<MCSUP>
      end
      
      % get.vtStimulusEndTimes - GETTER for 'vtStimulusEndTimes'
      function [vtStimulusEndTimes] = get.vtStimulusEndTimes(oStack)

         % - If the stimulus end times have been assigned, just return them
         if (oStack.bAssignedStimEndTimes)
            vtStimulusEndTimes = oStack.vtStimulusEndTimes;

         elseif (oStack.bAssignedStimStartTimes)
            % - Try to calculate vtStimulusEndTimes from start times and
            % duration
            vtStimDurations = [oStack.vtStimulusDurations; nan];
            cvnSequenceIDs = oStack.cvnSequenceIDs; %#ok<PROP>
            vnStimOrder = vertcat(cvnSequenceIDs{:}); %#ok<PROP>
            vnStimOrder(isnan(vnStimOrder)) = numel(vtStimDurations);
            vtStimulusEndTimes = oStack.vtStimulusStartTimes + vtStimDurations(vnStimOrder);

         else
            % - Try to calcluate end times from durations and order only
            [nul, vtStimulusEndTimes] = calc_stim_start_end_times(oStack);
         end
      end
      
      % set.mtStimulusUseTimes - SETTER for 'mtStimulusUseTimes'
      function set.mtStimulusUseTimes(oStack, mtStimulusUseTimes)
         if (~isnumeric(mtStimulusUseTimes) || (size(mtStimulusUseTimes, 1) ~= oStack.nNumStimuli) && ~isempty(mtStimulusUseTimes)) %#ok<MCSUP>
            error('FocusStack:InvalidArgument', ...
               '*** FocusStack/set.mtStimulusUseTimes: ''mtStimulusUseTimes'' must have [%d x 2] elements for this stack.', ...
               oStack.nNumStimuli); %#ok<MCSUP>
         end
         
         % - Assign mtStimulusUseTimes
         oStack.mtStimulusUseTimes = mtStimulusUseTimes;
      end
      
      % get.mtStimulusUseTimes - GETTER for 'mtStimulusUseTimes'
      function [mtStimulusUseTimes] = get.mtStimulusUseTimes(oStack)
         % - Has a set of use times been assigned?
         if (isempty(oStack.mtStimulusUseTimes))
            % - Does a set of durations exist?
            if (~isempty(oStack.vtStimulusDurations))
               % - Assume we should use the full assigned duration
               mtStimulusUseTimes = [zeros(oStack.nNumStimuli, 1) reshape(oStack.vtStimulusDurations, [], 1)];
               
            else
               % - Return an empty matrix -- there's not much we can assume
               mtStimulusUseTimes = [];
            end
            
         else
            % - Return the assigned use times
            mtStimulusUseTimes = oStack.mtStimulusUseTimes;
         end
      end
      

%% --- Save and load functions
      
      % saveobj - SAVE FUNCTION
      function oData = saveobj(oStack)
         % - Record useful values from the stack
         oData.bWritable = oStack.bWritable;
         oData.cstrFilenames = oStack.cstrFilenames;
         oData.mfFrameShifts = oStack.mfFrameShifts;
         oData.cmfBlankFrames = oStack.cmfBlankFrames;
         oData.mnAssignedBlankMeanFrames = oStack.mnAssignedBlankMeanFrames;
         oData.mnAssignedBlankStdFrames = oStack.mnAssignedBlankStdFrames;
         oData.vfBlackTrace = oStack.vfBlackTrace;
         oData.bSubtractBlack = oStack.bSubtractBlack;
         oData.bSubtractBlank = oStack.bSubtractBlank;
         oData.bConvertToDFF = oStack.bConvertToDFF;
         oData.vtStimulusDurations = oStack.vtStimulusDurations;
         oData.mtStimulusUseTimes = oStack.mtStimulusUseTimes;
         oData.tFrameDuration = oStack.tFrameDuration;
         
         if (oStack.bAssignedStimStartTimes)
            oData.vtStimulusStartTimes = oStack.vtStimulusStartTimes;
         end
         
         if (oStack.bAssignedStimEndTimes)
            oData.vtStimulusEndTimes = oStack.vtStimulusEndTimes;
         end
         
         if (oStack.bAssignedSequenceIDs)
            oData.cvnSequenceIDs = oStack.cvnSequenceIDs;
         end
      end
   end
   
   methods (Static = true)
      % loadobj - LOAD FUNCTION
      function oStack = loadobj(oData)
         % - Construct a new stack
         oStack = FocusStack(oData.cstrFilenames, oData.bWritable);
         
         % - Assign saved values
         if (~isempty(oData.mfFrameShifts))
            oStack.mfFrameShifts = oData.mfFrameShifts;
         end
         
         if (~isempty(oData.cmfBlankFrames))
            oStack.cmfBlankFrames = oData.cmfBlankFrames;
         end
         if (~isempty(oData.mnAssignedBlankMeanFrames))
            oStack.mnAssignedBlankMeanFrames = oData.mnAssignedBlankMeanFrames;
         end
         if (~isempty(oData.mnAssignedBlankStdFrames))
            oStack.mnAssignedBlankStdFrames = oData.mnAssignedBlankStdFrames;
         end
         if (~isempty(oData.vfBlackTrace))
            oStack.vfBlackTrace = oData.vfBlackTrace;
         end
         if (isfield(oData, 'cvnSequenceIDs') && ~isempty(oData.cvnSequenceIDs))
            oStack.cvnSequenceIDs = oData.cvnSequenceIDs;
         end
         if (isfield(oData, 'vtStimulusDurations') && ~isempty(oData.vtStimulusDurations))
            oStack.vtStimulusDurations = oData.vtStimulusDurations;
         end
         
         if (isfield(oData, 'vtStimulusStartTimes') && ~isempty(oData.vtStimulusStartTimes))
            oStack.vtStimulusStartTimes = oData.vtStimulusStartTimes;
         end
         
         if (isfield(oData, 'vtStimulusEndTimes') && ~isempty(oData.vtStimulusEndTimes))
            oStack.vtStimulusEndTimes = oData.vtStimulusEndTimes;
         end
         
         if (isfield(oData, 'bSubtractBlack'))
            oStack.bSubtractBlack = oData.bSubtractBlack;
         end
         
         if (isfield(oData, 'bSubtractBlank'))
            oStack.bSubtractBlank = oData.bSubtractBlank;
         end
         
         if (isfield(oData, 'bConvertToDFF'))
            oStack.bConvertToDFF = oData.bConvertToDFF;
         end
         
         if (isempty(oStack.tFrameDuration))
            oStack.tFrameDuration = oData.tFrameDuration;
         end
      end
   end
end

%% --- Helper functions

% calc_stim_start_end_time - FUNCTION Compute dynamic start and end times
function [vtStimulusStartTimes, vtStimulusEndTimes] = calc_stim_start_end_times(oStack)
            
   % - We need sequence IDs, stimulus durations and frame duration
   if (  isempty(oStack.vtStimulusDurations) || ...
         isempty(oStack.cvnSequenceIDs) || ...
         isempty(oStack.tFrameDuration))
      % - Just return an empty matrix;
      vtStimulusStartTimes = [];
      vtStimulusEndTimes = [];
      return;
   end

   % - Compute the start times
   cvnSequenceIDs = oStack.cvnSequenceIDs;
   vnStimOrder = vertcat(cvnSequenceIDs{:});
   nNumPresentations = numel(vnStimOrder);
   nFramesPerBlock = size(oStack, 3) / numel(oStack.cstrFilenames);
   
   tBlockDuration = (nFramesPerBlock * oStack.tFrameDuration);
   
   vtStimulusDurations = oStack.vtStimulusDurations; %#ok<PROP>
   nSkipIndex = numel(vtStimulusDurations)+1; %#ok<PROP>
   vtStimulusDurations(nSkipIndex) = nan; %#ok<PROP>
   vnStimOrder(isnan(vnStimOrder)) = nSkipIndex;
   
   vtStimulusStartTimes = [0; reshape(vtStimulusDurations(vnStimOrder), [], 1)]; %#ok<PROP>
   
   % - Replace NaNs with block skip durations
   while (any(isnan(vtStimulusStartTimes)))
      nSkipPresentation = find(isnan(vtStimulusStartTimes), 1, 'first');
      tCumulativeTime = sum(vtStimulusStartTimes(1:nSkipPresentation-1));
      tTimeLeftInBlock = tBlockDuration - mod(tCumulativeTime, tBlockDuration);
      vtStimulusStartTimes(nSkipPresentation) = tTimeLeftInBlock;
   end
   
   % - Compute start and end times
   vtStimulusEndTimes = cumsum(vtStimulusStartTimes); %#ok<PROP>
   vtStimulusEndTimes = vtStimulusEndTimes(2:end); %#ok<PROP>
   vtStimulusStartTimes = cumsum(vtStimulusStartTimes(1:end-1));
end

% --- END of FocusStack CLASSDEF ---
