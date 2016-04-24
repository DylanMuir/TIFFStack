function [cellArray] = CellFlatten(varargin)

% CellFlatten - FUNCTION Convert a list of items to a single level cell array
% $Id: CellFlatten.m 181 2005-02-27 19:07:59Z dylan $
%
% Usage: [cellArray] = CellFlatten(arg1, arg2, ...)
%
% CellFlatten will convert a list of arguments into a single-level cell array.
% If any argument is already a cell array, each cell will be concatenated to
% 'cellArray' in a list.  The result of this function is a single-dimensioned
% cell array containing a cell for each individual item passed to CellFlatten.
% The order of cell elements in the argument list is guaranteed to be
% preserved.
%
% This function is useful when dealing with variable-length argument lists,
% each item of which can also be a cell array of items.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 14th May, 2004

% -- Check arguments

if (nargin == 0)
   disp('*** CellFlatten: Do you want help?');
   help private/CellFlatten;
   return;
end


% -- Convert arguments

if (iscell(varargin{1}))
   cellArray = CellFlatten(varargin{1}{:});
else
   cellArray = {varargin{1}};
end

for (nIndexArg = 2:length(varargin))
   if (iscell(varargin{nIndexArg}))
      cellReturn = CellFlatten(varargin{nIndexArg}{:});
      cellArray = {cellArray{:} cellReturn{:}};
   else
      cellArray = {cellArray{:} varargin{nIndexArg}};
   end
end


% --- END of CellFlatten.m ---
