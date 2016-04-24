function [varargout] = subsref(oStack, sSubs)

% subsref - METHOD Overloaded subsref for focus stacks

if (isequal(sSubs(1).type, '.'))
   % - Handle referencing of high- and low-byte image stacks
   switch (sSubs(1).subs)
      case {'HighByteChannel'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', 2};
         end
         
         varargout{1} = ExtractFrames(oStack, sSubs(2));
   
      case {'LowByteChannel'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', 1};
         end
         
         varargout{1} = ExtractFrames(oStack, sSubs(2));

      case {'RawStack'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', ':'};
         end
         
         varargout{1} = ExtractFrames(oStack, sSubs(2));

      case {'AlignedStack'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', ':'};
         end
         
         varargout{1} = ExtractAlignedFrames(oStack, sSubs(2));
      
      case {'AlignedStackDouble'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', ':'};
         end
         
         varargout{1} = ExtractAlignedFramesDouble(oStack, sSubs(2));
         
         
      case {'SummedFrames'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', ':'};
         end
         
         varargout{1} = ExtractSummedFrames(oStack, sSubs(2), false);

      case {'SummedAlignedFrames'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', ':'};
         end
         
         varargout{1} = ExtractSummedFrames(oStack, sSubs(2), true);

      case {'BlankFrames'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':'};
         end

         [varargout{1:2}] = ExtractBlankFrames(oStack, sSubs(2));
         
      case {'MeanPixels'}
         if (numel(sSubs) == 1)
            sSubs(2).type = '()';
            sSubs(2).subs = {':', ':', ':', ':'};
         end
         
         varargout{1} = ExtractMeanPixels(oStack, sSubs(2), true);

      case {'FrameStimulusInfo'}
         [varargout{1:8}] = FrameStimulusInfo(oStack, sSubs(2).subs{1});
            
      otherwise
         % - Call the built-in subsref
         varargout{1} = builtin('subsref', oStack, sSubs);
   end
   
elseif (isequal(sSubs(1).type, '()'))
   % - Handle a direct reference
   if isAligned(oStack)
      varargout{1} = ExtractAlignedFrames(oStack, sSubs);
   else
      varargout{1} = ExtractFrames(oStack, sSubs);
   end
   
else
   % - This is an error
   error('FocusStack:InvalidReference', '*** FocusStack/subsref: Invalid reference.');
end

% --- END of subsref METHOD ---
