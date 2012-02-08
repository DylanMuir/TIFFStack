function CloseFiles(oStack)

% CloseFiles - METHOD Close the file handles for this stack

% This is taken care of by the memmapfile objects themselves.

for (hHandle = oStack.vhMemMapFileHandles(:)')
   switch (class(hHandle))
      case 'memmapfile'
         hHandle{1}.delete();
         
      case 'TIFFStack'
         hHandle.delete();
   end
end

if (~isempty(oStack.oAlignedFrameCache))
   delete(oStack.oAlignedFrameCache);
end

% --- END of CloseFiles METHOD ---
