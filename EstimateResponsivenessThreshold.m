function [fZThresh, fEquivAlpha] = EstimateResponsivenessThreshold(tfResponseZScores, sCellROIs, fDesiredAlpha)

% EstimateResponsivenessThreshold - FUNCTION Estimate a Z score threshold for a desired neuropil significance
%
% Usage: [fZThresh, fEquivAlpha] = EstimateResponsivenessThreshold(tfResponseZScores, sCellROIs, fDesiredAlpha)

if (nargin < 3)
   disp('*** EstimateResponsivenessThreshold: Incorrect usage.');
   help EstimateResponsivenessThreshold;
end

% - Make a non-ROI-mask
mbNeuropilMask = ~labelmatrix(sCellROIs);

% - Select non-ROI Z scores and sort
vfNeuropilZScores = tfResponseZScores(repmat(mbNeuropilMask, [1 1 size(tfResponseZScores, 3)]));
vfNeuropilZScores = vfNeuropilZScores(vfNeuropilZScores > 0);  % Accept only positive deviations
vfNeuropilZScores = sort(vfNeuropilZScores);

% - Find a threshold value that meets the desired alpha
fZThresh = vfNeuropilZScores(ceil(numel(vfNeuropilZScores)*(1-fDesiredAlpha)));
fEquivAlpha = 1-normcdf(fZThresh);

% --- END of EstimateResponsivenessThreshold.m ---

