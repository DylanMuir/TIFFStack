function opentif(strFilename)

% opentif - FUNCTION Drag-and-drop target for TIF stack files
%
% Usage: opentif(strFilename)
%
% This function acts as a file drop target for matlab to process TIF stack
% (.tif) files.  Just have it on your matlab path.

QuickProcessStack(strFilename);

% --- END of opentif.m ---
