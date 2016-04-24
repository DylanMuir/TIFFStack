function [mfCorr] = CorrelationImage(oStack, oActivity, vnFrameRange, nChannel, mfMeanAct)

% CorrelationImage - FUNCTION Measure the correlation of a pixel with the stack
%
% Usage: [mfCorr] = CorrelationImage(oStack, sSubs, [], [], mfAvgResponse)
%        [mtCorr] = CorrelationImage(oStack, mfActivityTrace, vnFrameRange, nChannel, mfMeanAct)
%                   CorrelationImage(oStack, cRefs, [], [], mfAvgResponse)
%
% If 'mfActivityTrace' is supplied as a matrix, then the correlation with
% several activity traces can be measured simultaneously.  In that case,
% 'mfActivityTrace' must be [NxT], where N is the number of activity traces.
% 'tfCorr' will be [XxYxN], where XxY are the dimensions of a single activity
% frame.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 25th March, 2011

% -- Check arguments

if (nargin < 2)
   disp('*** CorrelationImage: Incorrect usage');
   help CorrelationImage;
   return;
end


% -- Make sure we have an activity trace to correlate

subs = [];

% - Is this a subs structure?
if (isstruct(oActivity) && isfield(oActivity, 'subs'))
   subs = oActivity.subs;
   
elseif (iscell(oActivity))
   subs = oActivity;
end

% - Try to extract an activity trace
if (~isempty(subs))
   mfActivityTrace = permute(ExtractAlignedFrames(oStack, subs), [1 3 2]);
   nNumTraces = 1;

   % - Extract reference parameters
   [cvnFileRefs, cvnFullRefs] = GetFullFileRefs(oStack, subs);
   vnFrameRange = cvnFullRefs{1};
   nChannel = cvnFullRefs{2}(1);
   
else
   if (nargin < 4)
      disp('*** CorrelationImage: ''vnFrameRange'' and ''nChannel'' must be provided for a vector activity trace.');
      help CorrelationImage;
      return;
   end
   
   mfActivityTrace = oActivity;
   nNumTraces = size(mfActivityTrace, 1);
end


% - Extract the mean responses
vnStackSize = size(oStack);

if (~exist('mfMeanAct', 'var') || isempty(mfMeanAct))
   mfMeanAct = ExtractSummedFrames(oStack, {1:vnStackSize(1), 1:vnStackSize(2), vnFrameRange, nChannel}, true) ./ numel(vnFrameRange);
elseif (~isequal(size(mfMeanAct), vnStackSize(1:2)))
   error('FocusStack:InvalidArgument', ...
         '*** CorrelationImage: Error: ''mfMeanAct'' must be the same size as a stack frame.');
end

vfMeanTestAct = mean(mfActivityTrace, 2);
mfTestDeviation = double(mfActivityTrace) - repmat(vfMeanTestAct, 1, numel(vnFrameRange));


% -- Accumulate the correlation

tfNumerator = zeros([size(oStack, 1:2), nNumTraces]);
mfSumSqrDeviation = zeros(size(oStack, 1:2));

vnFrameOnes1 = ones(1, vnStackSize(1));
vnFrameOnes2 = ones(1, vnStackSize(2));
vnTracesOnes = ones([1 1 nNumTraces]);

mfBlockTestDeviation = ones([vnStackSize(1:2) nNumTraces]);
mfBlockThisDeviation = ones([vnStackSize(1:2) nNumTraces]);

for (nFrame = 1:numel(vnFrameRange))
   % - Get this frame
   mfThisFrame = ExtractAlignedFrames(oStack, {1:vnStackSize(1), 1:vnStackSize(2), vnFrameRange(nFrame), nChannel});
   mfThisDeviation = double(mfThisFrame) - mfMeanAct;
   
   % - Accumulate interaction term
   vfThisDeviation = permute(mfTestDeviation(:, nFrame), [2 3 1]);
   mfBlockTestDeviation(:, :, :) = vfThisDeviation(vnFrameOnes1, vnFrameOnes2, :);
   mfBlockThisDeviation(:, :, :) = mfThisDeviation(:, :, vnTracesOnes);
   tfNumerator = tfNumerator + mfBlockTestDeviation .* mfBlockThisDeviation;
   
   mfSumSqrDeviation = mfSumSqrDeviation + ...
      mfThisDeviation.^2;
end

% - Compute correlation
mfCorr = tfNumerator ./ (repmat(permute(sqrt(sum(mfTestDeviation.^2, 2)), [2 3 1]), vnStackSize(1:2)) .* repmat(sqrt(mfSumSqrDeviation), [1 1 nNumTraces]));

% --- END of CorrelationImage.m ---
