function [vfPRI, vbSignificantPRI] = ...
   PlaidVsComponentModel(  tfTrialResponses, ...
                           mfStimZScores, fSigZThresh, ...
                           mnStimuliDefinition, vnComponentMatchingPlaid)

% PlaidResponseIndex - FUNCTION Compute the plaid response index for a population
%
% Usage: [vfPRI, vbSignificantPRI] = ...
%           PlaidVsComponentModel(tfTrialResponses, ...
%                                 mfStimZScores, fSigZThresh, ...
%                                 mnStimuliDefinition, vnComponentMatchingPlaid)
%
% A PRI of 1 means the ROI is purely pattern-selective.  A PRI of -1 means the
% ROI is purely component-selective.

% Created: 6th October, 2011
% Author: Dylan Muir <muir@hifo.uzh.ch>


% -- Check inputs

nNumROIs = size(tfTrialResponses, 1);
nNumStimuli = size(tfTrialResponses, 2);
nNumTrials = size(tfTrialResponses, 3);
nNumComponentsInPlaid = size(mnStimuliDefinition, 2);

if (nargin < 5)
   disp('*** PlaidVsComponentModel: Incorrect usage');
end

if (~isequal(size(mfStimZScores), [nNumROI nNumStimuli]))
   disp('*** PlaidVsComponentModel: ''tfTrialResponses'' and ''mfStimZScores'' must have the same numbers of ROIs and stimuli.');
   return;
end

% -- Make sure everything is the right shape

vnComponentMatchingPlaid = reshape(vnComponentMatchingPlaid, [], 1);


% -- Tease apart stimulus definition

vbIsBlank = all(mnStimuliDefinition == 0, 2);
vbIsPlaid = all(mnStimuliDefinition > 0, 2);
vbIsComponent = ~vbIsPlaid & ~vbIsBlank;
vnComponentIndex = zeros(nNumStimuli, 1);
vnComponentIndex(vbIsComponent) = sum(mnStimuliDefinition(vbIsComponent, :), 2); % Because one entry is non-zero for each row

vnIndices = 1:nNumStimuli;
vnCompLookup = reshape(vnIndices(vbIsComponent), [], 1);


% -- Build the "component-selective" model for each ROI
% - "Component-selective" means that the response to a plaid is the summed
% responses to the two components of the plaid.

tnROIIndices = repmat((1:nNumROIs)', [1 nnz(vbIsPlaid) nNumComponentsInPlaid]);
tnComponentIndices = repmat(permute(vnCompLookup(mnStimuliDefinition(vbIsPlaid, :)), [3 1 2]), [nNumROIs 1 1]);

tfPlaidComponents = mfResponses(sub2ind(size(mfResponses), tnROIIndices, tnComponentIndices));
tfPlaidZScores = mfStimZScores(sub2ind(size(mfResponses), tnROIIndices, tnComponentIndices));

mfComponentSelectiveModel = sum(tfPlaidComponents, 3);



% -- Build the "plaid-selective" model for each ROI
% - "Plaid selective" means that the response to a plaid is equal to the
% response of the single component that moves in the same direction

mnROIIndices = repmat((1:nNumROIs)', [1 nnz(vbIsPlaid)]);
mnComponentIndices = repmat(reshape(vnCompLookup(vnComponentMatchingPlaid(vbIsPlaid)), 1, []), [nNumROIs 1]);
mfPlaidSelectiveModel = mfResponses(sub2ind(size(mfResponses), mnROIIndices, mnComponentIndices));


% -- Compute correlations between models and responses

vfCompSelCorr = diag(corr(mfResponses(:, vbIsPlaid)', mfComponentSelectiveModel'));
vfPlaidSelCorr = diag(corr(mfResponses(:, vbIsPlaid)', mfPlaidSelectiveModel'));
vfModelCorr = diag(corr(mfComponentSelectiveModel', mfPlaidSelectiveModel'));

vfCompSelPartCorr = (vfCompSelCorr - vfPlaidSelCorr .* vfModelCorr) ./ sqrt((1-vfPlaidSelCorr.^2) .* (1-vfModelCorr.^2));
vfPlaidSelPartCorr = (vfPlaidSelCorr - vfCompSelCorr .* vfModelCorr) ./ sqrt((1-vfCompSelCorr.^2) .* (1-vfModelCorr.^2));

% - Convert to z'-scores (Fisher's r-to-z' transformation)
%   z' scores have a mean of zero and a standard error of 1/sqrt(nDOF - 3)
vfCompSelZ = atanh(vfCompSelPartCorr);
vfPlaidSelZ = atanh(vfPlaidSelPartCorr);

% nDOF = (nNumTrials-1) * (nnz(vbIsPlaid) + nnz(vbIsComponent)) - 3;
% nDOF = (nnz(vbIsPlaid) + nnz(vbIsComponent)) - 3;
nDOF = (nNumTrials-1) * nnz(vbIsComponent) - 3;


% - Convert to P value estimates
vfPvsCPVal = 1-normcdf(vfPlaidSelZ-vfCompSelZ, 0, sqrt(2/nDOF));
vfCvsPPVal = 1-normcdf(vfCompSelZ-vfPlaidSelZ, 0, sqrt(2/nDOF));
vfPvs0PVal = 1-normcdf(vfPlaidSelZ, 0, 1 / sqrt(nDOF));
vfCvs0PVal = 1-normcdf(vfCompSelZ, 0, 1 / sqrt(nDOF));


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
