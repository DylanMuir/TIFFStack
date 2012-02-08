function nEndVal = end(oStack, nIndexDim, nNumIndices)

% end - METHOD Overloaded 'end' method

vnSize = size(oStack);

if (nIndexDim < nNumIndices)
   nEndVal = vnSize(nIndexDim);
   
else
   nEndVal = prod(vnSize(nIndexDim:end));
end

% --- END of end METHOD ---
