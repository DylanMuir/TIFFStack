function oStack = cat(nul, varargin) %#ok<INUSL,INUSD>

% cat - METHOD Overloaded concatenation method

nNumStacks = numel(varargin);

% - Handle no stacks
if (nNumStacks == 0)
   oStack = [];
   return;
end

% - Handle a single stack
if (nNumStacks == 1)
   oStack = varargin{1};
   return;
end


% -- Handle multiple stacks by creating a new stack

cstrFilenames = {};

for (nStack = 1:nNumStacks)
   cstrFilenames = [cstrFilenames(:); varargin{nStack}.cstrFilenames]; %#ok<NODEF>
end

oStack = FocusStack(cstrFilenames, varargin{1}.bWritable);


% --- END of cat METHOD ---
