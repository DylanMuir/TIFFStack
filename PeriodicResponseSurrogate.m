function [vfData, vtTime] = PeriodicResponseSurrogate( nNumPeriods, tSampleRate, ...
                                             f2ndHarmonicAmplitude, fTuningWidth, ...
                                             vfLFNoiseAmpStdCutoffHz, fShotNoiseAmp)

%% Generate the ideal response

% - Make a time trace vector
vtTime = linspace(0, nNumPeriods, tSampleRate * nNumPeriods);

% - Calculate the ideal response
vfResponseDistance = min(abs(bsxfun(@minus, vtTime', (0:nNumPeriods))), [], 2)';
vfResponseDistance2nd = min(abs(bsxfun(@minus, vtTime', (0:nNumPeriods-1)+0.5)), [], 2)';
vfResponse = normpdf(vfResponseDistance, 0, fTuningWidth/4) + ...
               f2ndHarmonicAmplitude * normpdf(vfResponseDistance2nd, 0, fTuningWidth/4);

%% Add noise
                              
vfShotNoise = fShotNoiseAmp + (randn(size(vtTime)) * sqrt(fShotNoiseAmp));

if (exist('vfLFNoiseAmpStdCutoffHz', 'var') && ~isempty(vfLFNoiseAmpStdCutoffHz))
   vfLFTrace = vfLFNoiseAmpStdCutoffHz(1) .* normrnd(0, vfLFNoiseAmpStdCutoffHz(2), size(vtTime));
   fSampFreq = 1/diff(vtTime(1:2));
   NFFT = 2^nextpow2(numel(vtTime));
   vfFreq = fSampFreq/2*linspace(0, 1, NFFT/2);
   vfFreq = [vfFreq vfFreq(end:-1:1)];
   vbFilter = vfFreq <= vfLFNoiseAmpStdCutoffHz(3);
   vfLFTrace = real(ifft(fft(vfLFTrace, NFFT) .* vbFilter, NFFT));
   vfLFTrace = vfLFTrace(1:numel(vtTime));
   vfLFTrace = vfLFTrace-min(vfLFTrace);
   vfLFTrace = vfLFNoiseAmpStdCutoffHz(1) .* vfLFTrace./max(vfLFTrace);

else
   vfLFTrace = zeros(size(vtTime));
end

vfData = vfResponse + vfShotNoise + vfLFTrace;

% --- END of PeriodicResponseSurrogate.m ---
