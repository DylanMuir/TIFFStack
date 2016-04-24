function opentif(strFilename)

% opentif - FUNCTION Drag-and-drop target for TIF stack files
%
% Usage: opentif(strFilename)
%
% This function acts as a file drop target for matlab to process TIF stack
% (.tif) files.  Just have it on your matlab path.

persistent cellFilenames tTimer;

% QuickProcessStack(strFilename);

% - Append this filename(s)
cellFilenames = [cellFilenames {strFilename}];

% - Either create or extend the timer
if isempty(tTimer)
   tTimer = timer('Name', 'FocusStack open execution timer', 'StartDelay', 0.5, 'TimerFcn', @(o,e)openfcs_makestack(strFilename));
   start(tTimer);
else
   stop(tTimer);
   start(tTimer);
end


   function openfcs_makestack(strFilename)
      
      [~, name] = fileparts(strFilename);
      
      strVarName = ['fs_' name];
      fs = FocusStack(cellFilenames);
      cellFilenames = {};
      assignin('base', strVarName, fs);
      
      fprintf('--- openfcs: Created FocusStack [''%s'']\n', strVarName);
      PlayStack(fs, 1:size(fs, 4));
      
      tTimer = [];
   end

end


% --- END of opentif.m ---
