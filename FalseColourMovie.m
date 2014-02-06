function FalseColourMovie(strFilename, fsStack, sRegions, vfDataRange, nChannel, mfAverage, mfColormap, strMovieType, bInvert, vbStimFrames, nStimMarkerSize)

% FalseColourMovie - FUNCTION Write a false-colour movie of stack activity
%
% Usage: FalseColourMovie(strFilename, fsStack, sRegions, vfDataRange <, nChannel, mfAverage, mfColormap, strMovieType, bInvert, vbStimFrames, nStimMarkerSize>)
%
% Generates a false-colour activity movie from a stack 'fsStack', showing
% only activity within the regions 'sRegions'.  TIFF and AVI movies can be
% written to disk, with the format specified by the extension of
% 'strFilename'.
%
% 'fsStack' can be a matlab tensor, or a FocusStack object.  In this second
% case, more sophisticated options are available.
%
% Activity will be mapped within the range
% ['vfDataRange(1)'..'vfDataRange(2)'].  Responses within this range will
% be mapped to a colour map, optionally given in 'mfColormap'.  By default
% the map 'jet' is used.
%
% The optional argument 'nChannel' defines which channel of 'fsStack' to
% use.
%
% The average stack image will be used as a background for the movie.
% This can optionally be provided in 'mfAverage', otherwise the average of
% 'fsStack' will be extracted.
%
% 'strMovieType' optionally defines how to manipulate responses for the
% movie.  By default ('raw'), the raw responses will be used.  The other
% options are 'dF/F', which computes the normalised change in flourescence
% (F - Fmean) / Fmean; and 'z-score', which computes dF/F and then
% normalises by the standard deviation.  In both these cases, 'fsStack'
% must be a FocusStack object, and the baseline and std. dev. frames must
% be assigned for each frame in the stack.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 6th February, 2013

% - Default arguments
DEF_mfColormap = jet(256);
DEF_strMovieType = 'raw';
DEF_nStimMarkerSize = 10;


%% -- Check arguments

if (nargin < 4)
   disp('*** FalseColourMovie: Incorrect usage');
   help FalseColourMovie;
end

if (~exist('mfColormap', 'var') || isempty(mfColormap))
   mfColormap = DEF_mfColormap;
end

if (~exist('strMovieType', 'var') || isempty(strMovieType))
   strMovieType = DEF_strMovieType;
end

if (~exist('nChannel', 'var') || isempty(nChannel))
   nChannel = 1;
end

if (~exist('bInvert', 'var') || isempty(bInvert))
   bInvert = false;
end


% - Get some information
nNumFrames = size(fsStack, 3);

% - Make sure the stack is un-normalised
if (isa(fsStack, 'FocusStack'))
   strNorm = fsStack.BlankNormalisation('none');
end

% - Do we need to obtain a stack average?
if (~exist('mfAverage', 'var') || isempty('mfAverage'))
   disp('--- FalseColourMovie: Extracting average response...');
   mfAverage = fsStack.SummedAlignedFrames(:, :, :, nChannel);
end

% - Normalise average and scale
mfAverage = mfAverage - min(mfAverage(:));
mfAverage = mfAverage ./ max(mfAverage(:)) * 0.5;

if (bInvert)
   mfAverage = 1 - mfAverage;
end

if (~exist('vbStimFrames', 'var'))
   vbStimFrames = false(1, nNumFrames);
end

if (~exist('nStimMarkerSize', 'var'))
   nStimMarkerSize = DEF_nStimMarkerSize;
end

%% -- What type of activity movie should we export?

switch(lower(strMovieType))
   case {'raw'}
      fhGetFrame = @(nFrame)(fsStack(:, :, nFrame, nChannel));
      
   case {'dff', 'df/f'}
      if (~isa(fsStack, 'FocusStack'))
         error('*** FalseColourMovie: DF/F movies are only supported for FocusStack objects.');
      end
      fhGetFrame = @GetDFFFrame;
      
   case {'z', 'zscore', 'z-score'}
      if (~isa(fsStack, 'FocusStack'))
         error('*** FalseColourMovie: Z-score movies are only supported for FocusStack objects.');
      end
      fhGetFrame = @GetZFrame;
end


%% -- What format should we write?

[~, ~, strFormat] = fileparts(strFilename);

% - Does the file exist?
if (exist(strFilename, 'file'))
   disp('*** FalseColourMovie: Warning: This file will be overwritten.');
   disp('       Press any key to continue.');
   pause;
   
   % - Delete the file
   delete(strFilename);
end

if (isempty(strFormat))
   strFormat = '.tif';
end

switch (lower(strFormat))
   case {'.tif', '.tiff'}
      fhWriteFrame = @(tfFrame)(imwrite(tfFrame, strFilename, 'WriteMode', 'append', 'Compression', 'none'));
      fhCloseFile = @()(true);
      
   case '.avi'
      hAVI = avifile(strFilename); %#ok<NASGU>
      
      fhWriteFrame = @(tfFrame)(eval('hAVI = addframe(hAVI, tfFrame);'));
      fhCloseFile = @()(eval('hAVI = close(hAVI);'));
      
   otherwise
      error('*** FalseColourMovie: Unsupported movie format.');
end


%% -- Export frames

% - Generate a base frame and a region mask
tfBaseIm = repmat(mfAverage, [1 1 3]);
mbRegionMask = labelmatrix(sRegions) > 0;

fprintf('--- FalseColourMovie: Exporting frames [%3.0f%%]', 0);

for (nFrame = 1:nNumFrames)
   % - Start with a blank (average) image
   tfThisIm = tfBaseIm;
   
   % - Get data for this frame
   mfThisFrame = double(fhGetFrame(nFrame));
   
   % - Mask off regions
   mfThisFrame(~mbRegionMask) = nan;
   
   % - Map to vfDataRange
   mfThisFrame = mfThisFrame - vfDataRange(1);
   mfThisFrame = mfThisFrame ./ vfDataRange(2);
   mfThisFrame(mfThisFrame < 0) = nan;
   mfThisFrame(mfThisFrame > 1) = 1;
   
   % - Convert to colormap indices
   mnIndices = round(mfThisFrame * (size(mfColormap, 1)-1))+1;
   
   % - Find which parts of the frame to replace
   mbThisMask = mfThisFrame > 0;
   
   % - Replace this image frame with colormap values
   tfThisIm(repmat(mbThisMask, [1 1 3])) = mfColormap(mnIndices(mbThisMask), :);
   
   % - Add a stimulus marker
   if (vbStimFrames(nFrame))
      tfThisIm(end-nStimMarkerSize:end, end-nStimMarkerSize:end, :) = 1;
   end
   
   % - Write this image frame
   fhWriteFrame(permute(tfThisIm, [2 1 3]));
   
   % - Show some progress
   fprintf(1, '\b\b\b\b\b\b[%3.0f%%]', nFrame / nNumFrames * 100);
end

fprintf(1, '\n');

% - Close the movie file
fhCloseFile();

   function mfZFrame = GetZFrame(nFrame)
      % - Get raw frame
      mfRawFrame = double(fsStack(:, :, nFrame, nChannel));
      
      % - Convert to DF/F and divide by std. dev.
      [mfBaseline, mfStdDev] = fsStack.BlankFrames(:, :, nFrame);
      mfZFrame = (mfRawFrame ./ mfBaseline - 1) ./ mfStdDev;
   end

   function mfDFFFrame = GetDFFFrame(nFrame)
      % - Get raw frame
      mfRawFrame = double(fsStack(:, :, nFrame, nChannel));

      % - Convert to DF/F
      mfBaseline = fsStack.BlankFrames(:, :, nFrame);
      mfDFFFrame = mfRawFrame ./ mfBaseline - 1;
   end

% -  Restore normalisation, if necessary
if (isa(fsStack, 'FocusStack'))
   fsStack.BlankNormalisation(strNorm);
end

end

% --- END of FalseColourMovie.m ---




