% SCRIPT Generate surrogate calcium data

%% Parameters

vnImageSize = [40 40];
nNumROIs = 4;
nROIWidth = 6;

fSpontaneousRate = 0.1; % Hz

fSampleRate = 50;        % Hz
nNumTrials = 5;
tBlankTime = 5;         % sec
tSignalTime = 5;        % sec
mfResponses = [1 1  .5 1 0;
                0 .5 .25 1 2;
                1 1  0 0 0;
                1 1  0.5 0 0.5];   % Hz

% - Calcium signal parameters
tTransientTau = 2;      % sec
fROIAmplitudeStd = 1;
fShotNoiseStd = .4;
vfStaticNoiseMeanStd = [2 0.05];
fNeuropilContamination = 0.1;
vfLFNoiseAmpStdCutoffHz = [.1 1 1/10];
vfFBaselineMeanStd = [1 .1];


%% Make some ROIs

mnROILocs = rand(nNumROIs, 2) .* repmat(vnImageSize-nROIWidth, nNumROIs, 1) + nROIWidth/2;
vnX = 0:(vnImageSize(2)-1);
vnY = 0:(vnImageSize(2)-1);
[mnY, mnX] = meshgrid(vnX, vnY);

sROIs = bwconncomp(true(vnImageSize));
sROIs.NumObjects = nNumROIs;

for (nROIIndex = nNumROIs:-1:1)
   mnD = sqrt((mnX-mnROILocs(nROIIndex, 2)).^2 + (mnY-mnROILocs(nROIIndex, 1)).^2);
   sROIs.PixelIdxList{nROIIndex} = find(mnD <= nROIWidth/2);
end


%% Make a signal trace

tTrialDuration = (tBlankTime + tSignalTime) * size(mfResponses, 2);
vtTimeTrace = 0:(1/fSampleRate):tTrialDuration;

mfInstFiringRate = zeros(size(mfResponses, 1), numel(vtTimeTrace)) + fSpontaneousRate;
vnSignalTrace = zeros(1, numel(vtTimeTrace));

for (nRespIndex = 1:numel(vfResponses))
   tFirstTime = (nRespIndex-1)*(tBlankTime+tSignalTime);
   vbInSignal = (vtTimeTrace >= (tFirstTime+tBlankTime)) & (vtTimeTrace <= (tFirstTime+tBlankTime+tSignalTime));
   mfInstFiringRate(:, vbInSignal) = repmat(mfResponses(:, nRespIndex), 1, nnz(vbInSignal));
   vnSignalTrace(vbInSignal) = nRespIndex;
end


%% Generate calcium signal

[cvtSpikeTimes, ctfResponseTrace, cvfROIGain] = ...
   GenerateSurrogateCalciumSignal(  vtTimeTrace, mfInstFiringRate, nNumTrials, tTransientTau, sROIs, ...
                                    vfFBaselineMeanStd, fROIAmplitudeStd, fShotNoiseStd, vfStaticNoiseMeanStd, ...
                                    fNeuropilContamination, vfLFNoiseAmpStdCutoffHz);


%% Play a movie of each trial

for (nTrial = 1:nNumTrials)
   PlayStack(double(ctfResponseTrace{nTrial}));
end

%% Export tiffs

w = warning('off', 'MATLAB:DELETE:FileNotFound');

% - Export trial responses
for (nTrial = 1:nNumTrials)
   % - Remove an existing file
   delete(['channel-1-trial-' num2str(nTrial) '.tif']);
   
   for (nFrame = 1:numel(vtTimeTrace))
      imwrite(ctfResponseTrace{nTrial}(:, :, nFrame), ['channel-1-trial-' num2str(nTrial) '.tif'], 'tif', 'WriteMode', 'append');
   end
end

% - Export mask image
imwrite(labelmatrix(sROIs)>0, 'region-mask.tif');

% - Export mat data
save 'response-traces.mat' vnSignalTrace vtTimeTrace cvtSpikeTimes cvfROIGain;
save 'all-data.mat';

warning(w);
