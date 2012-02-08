function [vfCaTemplate] = CalciumTransientTemplate(tSamplingRate, tDuration, vtTaus)

% CalciumTransientTemplate - FUNCTION Generate a template calcium transient
%
% Usage: [vfCaTemplate] = CalciumTransientTemplate(tSamplingRate, tDuration, vtTaus)

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 23rd December, 2010

% -- Defaults

DEF_vtTaus = 5;


% -- Check arguemtns

if (nargin < 1)
   disp('*** CalicumTransientTemplate: Incorrect usage');
   help CalciumTransientTemplate;
   return;
end

if (~exist('vtTaus', 'var') || isempty(vtTaus))
   vtTaus = DEF_vtTaus;
end

if (~exist('tDuration', 'var') || isempty(tDuration))
   tDuration = 5 * max(vtTaus(:));
end

nNumExps = numel(vtTaus);


% -- Compute exponentials and sum

% - Compute time trace
tTime = 0:(1/tSamplingRate):tDuration;

vfCaTemplate = zeros(1, numel(tTime));

for (nExp = 1:nNumExps) %#ok<FORPF>
   vfCaTemplate = vfCaTemplate + exp(-tTime / vtTaus(nExp));
end

% - Normalise so that max response is 1
vfCaTemplate = vfCaTemplate ./ max(vfCaTemplate);

% --- END of CalciumTransientTemplate.m ---
