function [vfPRI, vbSignificantPRI] = ...
   PlaidResponseIndex(  mfResponses, ...
                        mfStimZScores, fSigZThresh, ...
                        vnComponentStimIDs, mnMatchingPlaidComparisonIDs)

% PlaidResponseIndex - FUNCTION Compute the plaid response index for a population
%
% Usage: [vfPRI, vbSignificantPRI] = ...
%           PlaidResponseIndex(mfResponses, ...
%                              mfStimZScores, fSigZThresh, ...
%                              vnComponentStimIDs, mnMatchingPlaidComparisonIDs)
%
% A PRI of 1 means the ROI is purely pattern-selective.  A PRI of -1 means the
% ROI is purely component-selective.

% Created: 6th October, 2011
% Author: Dylan Muir <muir@hifo.uzh.ch>


% -- Check inputs

nNumROIs = size(mfResponses, 1);
nNumStimuli = size(mfResponses, 2);

if (nargin < 5)
   disp('*** PlaidResponseIndex: Incorrect usage');
end

if (any(vnComponentStimIDs > nNumStimuli) || any(mnMatchingPlaidComparisonIDs(:) > nNumStimuli) || ...
    any(vnComponentStimIDs < 1) || any(mnMatchingPlaidComparisonIDs(:) < 1))
   disp('*** PlaidResponseIndex: The indices in ''vnComponentStimIDs'' and ''mnMatchingPlaidComparisonIDs'' must index into ''mfResponses''');
   return;
end

if (~isequal(size(mfResponses), size(mfStimZScores)))
   disp('*** PlaidResponseIndex: ''mfResponses'' and ''mfStimZScores'' must have the same size.');
   return;
end


% -- Compute the PRI for each ROI

% - Find max component response

[nul, vnMaxComponentID] = max(mfResponses(:, vnComponentStimIDs), [], 2);

% - Find corresponding plaid responses for analysis

mnStimCheckIndices = mnMatchingPlaidComparisonIDs(vnMaxComponentID, :);

% - Extract corresponding plaid responses
mnROIIndices = repmat(reshape(1:nNumROIs, [], 1), 1, 3);
vnPlaidIndices = sub2ind(size(mfResponses), mnROIIndices, mnStimCheckIndices);
mfPlaidResponses = mfResponses(vnPlaidIndices);

% - Compute PRI
% PRI = (1-(2+3)/(1+(2+3))
vfPRI = (mfPlaidResponses(:, 1) - sum(mfPlaidResponses(:, 2:3), 2)) ./ ...
         (mfPlaidResponses(:, 1) + sum(mfPlaidResponses(:, 2:3), 2));

% - Classify responses
vbPatternSel = vfPRI > 0;
vbComponentSel = vfPRI < 0;
      
% - Check for significant responses
vnComponentIndices = sub2ind(size(mfStimZScores), 1:nNumROIs, vnMaxComponentID');
vbSignificantMaxComponent = (mfStimZScores(vnComponentIndices) > fSigZThresh)';

vnPrincipalPlaidIndices = sub2ind(size(mfStimZScores), 1:nNumROIs, mnStimCheckIndices(:, 1)');
vbSignificantPrincipalPlaid = (mfStimZScores(vnPrincipalPlaidIndices) > fSigZThresh)';

vnOrthPlaid1Indices = sub2ind(size(mfStimZScores), 1:nNumROIs, mnStimCheckIndices(:, 2)');
vbSignificantOrth1Plaid = (mfStimZScores(vnOrthPlaid1Indices) > fSigZThresh)';

vnOrthPlaid2Indices = sub2ind(size(mfStimZScores), 1:nNumROIs, mnStimCheckIndices(:, 2)');
vbSignificantOrth2Plaid = (mfStimZScores(vnOrthPlaid2Indices) > fSigZThresh)';

% - Check for significance differently for pattern and component ROIs
vbSignificantPRI = false(nNumROIs, 1);
vbSignificantPRI(vbPatternSel) = vbSignificantMaxComponent(vbPatternSel) & vbSignificantPrincipalPlaid(vbPatternSel);
vbSignificantPRI(vbComponentSel) = vbSignificantMaxComponent(vbComponentSel) & vbSignificantOrth1Plaid(vbComponentSel) & vbSignificantOrth2Plaid(vbComponentSel);

% --- END of PlaidResponseIndex.m ---
