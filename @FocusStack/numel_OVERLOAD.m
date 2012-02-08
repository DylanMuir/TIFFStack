function nNumElems = numel_OVERLOAD(oStack, varargin)

% numel - METHOD Overloaded numel method

if (isempty(varargin))
   nNumElems = prod(size(oStack)); %#ok<PSIZE>

else
   nNumElems = 1;
   
end


% --- END of numel METHOD ---

