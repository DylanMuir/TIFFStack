function [vnSize] = size(oStack, nDim)

% size - METHOD Overloaded 'size' method

% - Determine total stack size
vnSize = [oStack.vnFrameSize sum([oStack.vnNumFrames]) oStack.nNumChannels];

if (exist('nDim', 'var') && isnumeric(nDim))
   if (nDim < 1)
      % - nDim is not within range
      error('MATLAB:getdimarg:dimensionMustBePositiveInteger', 'Dimension argument must be a positive integer scalar within indexing range.');
   
   elseif (nDim > 4)
      % - All trailing dimensions are of size 1
      vnSize = 1;
   
   else
      % - Return the requested index
      vnSize = vnSize(nDim);
   end
end


% --- END of size METHOD ---
