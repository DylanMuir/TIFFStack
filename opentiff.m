function opentiff(strFilename)

% opentiff - FUNCTION Drag-and-drop target for TIFF stack files
%
% Usage: opentiff(strFilename)
%
% This function acts as a file drop target for matlab to process TIFF stack
% (.tiff) files.  Just have it on your matlab path.

QuickProcessStack(strFilename);

% --- END of opentiff.m ---
