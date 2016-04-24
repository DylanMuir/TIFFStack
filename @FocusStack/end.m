function nEndVal = end(oStack, nIndexDim, nNumIndices)

% end - METHOD Overloaded 'end' method

vnSize = size(oStack);

switch (nNumIndices)
   case 1
      nEndVal = prod(vnSize);
      
   case 2
      error('FocusStack:InvalidReference', ...
            '*** FocusStack/End: Two-index referencing is not supported.');
         
   case 3
      switch (nIndexDim)
         case 1
            nEndVal = prod(vnSize(1:2));
            
         case 2
            nEndVal = vnSize(3);
            
         case 3
            nEndVal = vnSize(4);
            
         otherwise
            nEndVal = prod(vnSize((nIndexDim+1):end));
      end
      
   otherwise
      if (nIndexDim < nNumIndices)
         nEndVal = vnSize(nIndexDim);
         
      else
         nEndVal = prod(vnSize(nIndexDim:end));
      end
end


% --- END of end METHOD ---
