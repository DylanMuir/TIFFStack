function [vbResponsive] = IsResponsive(tfTrialResponsesDFF, mfStimZScore, fResponseThreshold, fZScoreThreshold)

nNumTrials = size(tfTrialResponsesDFF, 3);

vbZThresh = any(mfStimZScore > fZScoreThreshold, 2);

tbRespThresh = tfTrialResponsesDFF > fResponseThreshold;
mbReliableStim = sum(tbRespThresh, 3) > nNumTrials/2;

vbResponsive = vbZThresh & any(mbReliableStim, 2);
