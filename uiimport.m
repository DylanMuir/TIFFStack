function uiimport(varargin)

% uiimport - FUNCTION Drag-and-drop target for TIF stack files
%
% Usage: uiimport(varargin)
%
% This function acts as a file drop target for matlab to process TIF stack
% (.tif) files.  Just have it on your matlab path.

[nul, nul, strExt] = fileparts(varargin{1});

switch lower(strExt)
   case {'.tif', '.tiff'}
      QuickProcessStack(varargin{1});
      
   otherwise
      strPWD = cd(fullfile(matlabroot, 'toolbox', 'matlab', 'codetools'));
      uiimport(varargin{:});
      cd(strPWD);
end



% --- END of opentif.m ---
