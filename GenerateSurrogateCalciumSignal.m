function [cvtSpikeTimes, ctfResponseTrace, cvfROIGain, vfROIBaseline] = ...
   GenerateSurrogateCalciumSignal(  vtTimeTrace, mfInstSpikeRate, nNumTrials, vtTransientTau, sROI, ...
                                    vfFBaselineMeanStd, fROIAmplitudeStd, fShotNoiseStd, vfStaticNoiseMeanStd, ...
                                    fNeuropilContamination, vfLFNoiseAmpStdCutoffHz)

nNumROIs = sROI.NumObjects;

% - Duplicate a spike rate function for all ROIs, if only one was provided
if (size(mfInstSpikeRate, 1) == 1)
   mfInstSpikeRate = repmat(mfInstSpikeRate, nNumROIs, 1);
end

if (numel(vtTransientTau) == 1)
   vtTransientTau = repmat(vtTransientTau, nNumROIs, 1);
end



%% - Set up global parameters

vnY = 1:sROI.ImageSize(1);
vnX = 1:sROI.ImageSize(2);
[mnY, mnX] = meshgrid(vnY, vnX);

tTrainLength = max(vtTimeTrace) - min(vtTimeTrace);

   
%% - Produce a set of gain values for each ROI

for (nROIIndex = nNumROIs:-1:1)
   % - Find centre of ROI and jitter
   [vnY, vnX] = ind2sub(sROI.ImageSize, sROI.PixelIdxList{nROIIndex});
   vnROICentre = [mean(vnY) mean(vnX)];
   vnROICentre = vnROICentre + (rand(1, 2)-0.5) * 2;
   vfD = sqrt((vnY-vnROICentre(1)).^2 + (vnX-vnROICentre(2)).^2);
   
   % - Generate a Gaussian gain field for this ROI
   fROIWidth = sqrt(max((vnY - vnROICentre(1)).^2) + max((vnX - vnROICentre(2)).^2))*2;
   fROIAmp = normrnd(1, fROIAmplitudeStd);
   vfGains = fROIAmp * normpdf(vfD, 0, fROIWidth/4) ./ normpdf(0, 0, fROIWidth/4);
   
   % - Store gain for this ROI
   cvfROIGain{nROIIndex} = vfGains; %#ok<AGROW>
   
   % - Make a baseline for this ROI
   vfROIBaseline(nROIIndex) = normrnd(vfFBaselineMeanStd(1), vfFBaselineMeanStd(2)); %#ok<AGROW>
end


%% -- Produce a global static noise field

vfStaticNoise = normrnd(vfStaticNoiseMeanStd(1), vfStaticNoiseMeanStd(2), 1, prod(sROI.ImageSize));


%% -- Produce a set of trials with different spike instantiations and different noise
for (nTrialIndex = nNumTrials:-1:1)
   mfTrialCalciumSignalTrace = zeros([prod(sROI.ImageSize) numel(vtTimeTrace)]);
   
   for (nROIIndex = nNumROIs:-1:1)
      % - Generate a poisson train at the maximum spike rate, by drawing ISIs
      fMaxFreq = max(mfInstSpikeRate(nROIIndex, :));
      nISIsToGenerate = max(fix(fMaxFreq * tTrainLength * 1.3), 1) + 1;

      tCurrentTime = min(vtTimeTrace);
      vISIs = [];
      while (tCurrentTime < max(vtTimeTrace))
         vTheseISIs = exprnd(1/fMaxFreq, 1, nISIsToGenerate);
         vISIs = [vISIs vTheseISIs]; %#ok<AGROW>
         tCurrentTime = tCurrentTime + sum(vISIs);
      end
      
      % - Convert to spike times
      vtSpikeTimes = cumsum(vISIs);
      
      % - Pick out the spike times which fall inside the vTimeTrace window
      vbSpikesInTime = vtSpikeTimes <= max(vtTimeTrace);
      vtSpikeTimes = vtSpikeTimes(vbSpikesInTime);
      
      % - Work out instatnaeous spiking probability
      vfInstSpikeProb = interp1(vtTimeTrace, mfInstSpikeRate(nROIIndex, :), vtSpikeTimes, 'nearest') ./ fMaxFreq;
      
      % - Thin the spike train by only keeping spikes with high spiking probability
      vbKeepSpike = vfInstSpikeProb >= rand(1, numel(vfInstSpikeProb));
      vtSpikeTimes = vtSpikeTimes(vbKeepSpike);
      
      % - Record this spike train for this ROI
      cvtSpikeTimes{nROIIndex, nTrialIndex} = vtSpikeTimes; %#ok<AGROW>
      
      % - Convolve with transient function and accumulate
      vfROISignal = zeros(size(vtTimeTrace));
      for (nSpike = 1:numel(vtSpikeTimes))
         vtDiffTime = vtTimeTrace - vtSpikeTimes(nSpike);
         vtDiffTime(vtDiffTime<0) = inf;
         vfThisTransient = exp(-vtDiffTime) / vtTransientTau(nROIIndex);
         vfROISignal = vfROISignal + vfThisTransient;
      end
      
      % - Record this signal trace for this ROI
      cvfSignalTrace{nROIIndex, nTrialIndex} = vfROISignal; %#ok<AGROW>

      % - Accumulate into calcium signal trace
      mfTrialCalciumSignalTrace(sROI.PixelIdxList{nROIIndex}, :) = mfTrialCalciumSignalTrace(sROI.PixelIdxList{nROIIndex}, :) + ...
         repmat(cvfROIGain{nROIIndex}, 1, numel(vtTimeTrace)) .* repmat(vfROISignal, numel(sROI.PixelIdxList{nROIIndex}), 1) + vfROIBaseline(nROIIndex);
      
      %% - Add neuropil contamination for this ROI
      mfTrialCalciumSignalTrace = mfTrialCalciumSignalTrace + ...
         fNeuropilContamination .* repmat(vfROISignal, prod(sROI.ImageSize), 1);
   end
   
   %% - Generate pixel shot noise
   if (exist('fShotNoiseStd', 'var') && ~isempty(fShotNoiseStd))
      mfTrialCalciumSignalTrace = mfTrialCalciumSignalTrace + normrnd(0, fShotNoiseStd, size(mfTrialCalciumSignalTrace));
   end
   
   %% - Generate global low-frequency noise
   if (exist('vfLFNoiseAmpStdCutoffHz', 'var') && ~isempty(vfLFNoiseAmpStdCutoffHz))
      vfLFTrace = vfLFNoiseAmpStdCutoffHz(1) .* normrnd(0, vfLFNoiseAmpStdCutoffHz(2), size(vtTimeTrace));
      fSampFreq = 1/diff(vtTimeTrace(1:2));
      NFFT = 2^nextpow2(numel(vtTimeTrace));
      vfFreq = fSampFreq/2*linspace(0, 1, NFFT/2);
      vfFreq = [vfFreq vfFreq(end:-1:1)];
      vbFilter = vfFreq <= vfLFNoiseAmpStdCutoffHz(3);
      vfLFTrace = real(ifft(fft(vfLFTrace, NFFT) .* vbFilter, NFFT));
      vfLFTrace = vfLFTrace(1:numel(vtTimeTrace));
      vfLFTrace = vfLFTrace-min(vfLFTrace);
      vfLFTrace = vfLFNoiseAmpStdCutoffHz(1) .* vfLFTrace./max(vfLFTrace);
      
      % - Apply global LF noise trace
      mfTrialCalciumSignalTrace = mfTrialCalciumSignalTrace + repmat(vfLFTrace, prod(sROI.ImageSize), 1);
   end
   
   %% - Add static noise
   mfTrialCalciumSignalTrace = mfTrialCalciumSignalTrace + repmat(vfStaticNoise', 1, numel(vtTimeTrace));

   
   %% - Clip trace and discretise
   
   mfTrialCalciumSignalTrace(mfTrialCalciumSignalTrace<0)=0;
   mfTrialCalciumSignalTrace = uint8(round(mfTrialCalciumSignalTrace ./ max(mfTrialCalciumSignalTrace(:)) * 255));
   
   %% - Store this response trace
   ctfResponseTrace{nTrialIndex} = reshape(mfTrialCalciumSignalTrace, sROI.ImageSize(1), sROI.ImageSize(2), numel(vtTimeTrace)); %#ok<AGROW>
end

% --- END of GenerateSurrogateCalciumSignal.m ---
