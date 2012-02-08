function [vfPreferredAngle, vfMaxResponse, vfSelectivityIndex, vfVMKappa] = ...
   AnalyseOriDirSelectivity(fsData, sAnalysis, vfStimAngles) %#ok<INUSL>

% AnalyseOriDirSelectivity - FUNCTION Analyse orientation and direction selectivity of a population
%
% Usage: [vfPreferredAngle, vfMaxResponse, vfSelectivityIndex, vfVMKappa] = ...
%    AnalyseOriDirSelectivity(fsData, sAnalysis, vfStimAngles)

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 19th January, 2011

% -- Check arguments

if (nargin < 3)
   disp('--- AnalyseOriDirSelectivity: Incorrect usage');
   help AnalyseOriDirSelectivity;
   return;
end


% -- Analyse each region

vfInterpAngles = linspace(0, 2*pi, 200);

for (nRegion = sAnalysis.sRespRegions.NumObjects:-1:1)
   % - Ignore blank stimulus (stimulus 1)
   vfRegionResponse = sAnalysis.mfStimMeans(nRegion, 2:end);
   vfRegionResponse = vfRegionResponse ./ sum(vfRegionResponse);
   [vfPreferredAngle(nRegion), nul, vfVMKappa(nRegion)] = circ_mean_var(vfStimAngles, vfRegionResponse); %#ok<AGROW>
   
   % - Make a plot of this response
   figure;
   plot(vfStimAngles, vfRegionResponse, 'k-', 'LineWidth', 2);
   hold on;
   plot(vfInterpAngles, vonmisespdf(vfInterpAngles, vfPreferredAngle(nRegion), vfVMKappa(nRegion)), 'b:');
end

vfSelectivityIndex = [];
vfMaxResponse = [];

% --- END of AnalyseOriDirSelectivity.m ---
