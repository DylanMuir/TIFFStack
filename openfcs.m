function openfcs(strFilename)

% openfcs - FUNCTION Drag-and-drop target for FocusStack files
%
% Usage: openfcs(strFilename)
%
% This function acts as a file drop target for matlab to process Focus
% (.fcs) files.  Just have it on your matlab path.

% QuickProcessStack(strFilename);

[~, name] = fileparts(strFilename);

strVarName = ['fs_' name];
fs = FocusStack(strFilename);
assignin('base', strVarName, fs);

fprintf('--- openfcs: Created FocusStack [''%s'']\n', strVarName);
PlayStack(fs, 1:size(fs, 4));

% --- END of openfcs.m ---
