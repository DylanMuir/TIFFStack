function [fPopMean, fPopStd] = pop_mean_std(tfMeans, tfStds, tfWeights, nDim, bPopulation)

% pop_mean_std - FUNCTION Compute the population mean and std. dev. from sample means and std. devs
%
% Usage: [fPopMean, fPopStd] = pop_mean_std(tfMeans, tfStds, tfWeights <, nDim, bPopulation>)
%
% 'tfMeans' and 'tfStds' are vectors of mean and standard deviation estimates,
% respectively.
%
% 'fPopMean' and 'fPopStd' will be estimates for the aggregated mean and
% standard deviations for the combined population.
%
% 'tfWeights' allows sample sizes, corresponding to each element in 'vfMeans'
% and 'vfStds' to be provided.  For example, if 'vfMeans(1)' and 'vfStds(1)'
% were measured from 100 samples, but 'vfMeans(2)' and 'vfStds(2)' were measured
% from only 50 samples, then 'vfWeights' should be [100 50].  This corrects for
% different sample sizes.
%
% By default, sample-based statistics are used -- it is assumed that the entire
% population was NOT sampled.  To use population-based statistics, set the
% optional argument 'bPopulation' to true.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 22nd February, 2011

% -- Check arguments

if (nargin < 3)
   disp('*** pop_mean_std: Incorrect usage');
   help pop_mean_std;
   return;
end

if (~exist('nDim', 'var') || isempty(nDim))
   nDim = find(size(tfMeans) > 1, 1, 'first');
end

nNumSamples = size(tfMeans, nDim);
vnDataSize = size(tfMeans);

if (~isequal(size(tfStds), vnDataSize))
   disp('*** pop_mean_std: ''tfStds'' must have the same number of elements as ''tfMeans''');
   return;
end

if (~isequal(size(tfWeights), vnDataSize))
   disp('*** pop_mean_std: ''tfWeights'' must have the same number of elements as ''tfMeans'' and ''tfStds''');
   return;
end

if (~exist('bPopulation', 'var') || isempty(bPopulation))
   bPopulation = false;
end

% - Determine sample or population-based statistics
if (bPopulation)
   nSubtract = 0;
else
   nSubtract = 1;
end


% -- Combine samples

fPopMean = sum(tfWeights .* tfMeans, nDim) ./ sum(tfWeights, nDim);
fPopStd = sqrt(1 ./ (sum(tfWeights, nDim)-nSubtract) .* ...
               (sum((tfWeights-nSubtract) .* (tfStds.^2) + tfWeights .* (tfMeans.^2), nDim) - (sum(tfWeights, nDim) .* mean(tfMeans, nDim).^2)));

% --- END of pop_mean_std.m ---
