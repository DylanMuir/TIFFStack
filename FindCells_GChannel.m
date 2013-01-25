function [sCellRegions] = FindCells_GChannel(fsData, fLFRadiusUM, fHFRadiusUM, fCellRadiusUM)

% FindCells_GRChannels - FUNCTION Try to locate neuron somata using only the green channel
%
% Usage: [sCellRegions] = FindCells_GChannel(fsData <, fLFRadiusUM, fHFRadiusUM, fCellRadiusUM>)
%
% 'fsData' is a FocusStack object.  The average images (taken from stimulus
% response frames only, if the stimulation times are available in the
% stack) will be used. Band-pass filtering will be used to locate neurons.
%
% 'sCellRegions' will be a region structure, as returned by matlab function
% bwconncomp, defining the circular regions corresponding to each cell.
%
% The optional parameters 'fLFRadiusUM', 'fHFRadiusUM' and 'fCellRadiusUM'
% control the filtering process.  By default, they adopt reasonable values for
% sulfarhodamine plus injected OGB labelling in mouse visual cortex.  All values
% are in micrometres, and rely on the spatial calibration of the stack.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: November, 2011

% -- Default parameters

DEF_fCellRadiusUM = 6;
DEF_fHFRadiusUM = 4;
DEF_fLFRadiusUM = 25;


% -- Check arguments

if (~exist('fLFRadiusUM', 'var') || isempty(fLFRadiusUM))
   fLFRadiusUM = DEF_fLFRadiusUM;
end

if (~exist('fHFRadiusUM', 'var') || isempty(fHFRadiusUM))
   fHFRadiusUM = DEF_fHFRadiusUM;
end

if (~exist('fCellRadiusUM', 'var') || isempty(fCellRadiusUM))
   fCellRadiusUM = DEF_fCellRadiusUM;
end

%% -- Calibrate filtering parameters

if (isempty(fsData.fPixelsPerUM))
   warning('FocusStack:Uncalibrated', ...
      '--- FocusStack/FindCells_GRChannels: The FocusStack object does not contain calibration information.  Assuming 1 pixel = 1 um.');
   fPixelsPerUM = 1;
else
   fPixelsPerUM = fsData.fPixelsPerUM;
end

fLFRadiusPix = fLFRadiusUM * fPixelsPerUM;
fHFRadiusPix = fHFRadiusUM * fPixelsPerUM;
fCellRadiusPix = fCellRadiusUM * fPixelsPerUM;


%% -- Filter average frames

% - Turn of DFF conversion
strNormalisation = fsData.BlankNormalisation('none');

% - Get frames for stimulus responses for green channel
[  vtGlobalTime, ...
   vnBlockIndex, vnFrameInBlock, vtTimeInBlock, ...
   vnStimulusSeqID, vtTimeInStimPresentation, ...
   vnPresentationIndex, vbUseFrame] = ...
      FrameStimulusInfo(fsData, 1:size(fsData, 3));

% - If no stimuli were defined, use all frames
vnStackSize = size(fsData);
if (isempty(vbUseFrame))
   vbUseFrame = true(1, vnStackSize(3));
end

% - Trun off alignment warning
wOld = warning('off', 'FocusStack:UnalignedStack');

% - Get average green channel response
mfG = double(fsData.SummedAlignedFrames(:, :, vbUseFrame, 1)) ./ vnStackSize(3);

% - Construct filters, and filter average frames
hGLF = fspecial('gaussian', ceil(fLFRadiusPix*5), fLFRadiusPix);
hGHF = fspecial('gaussian', ceil(fLFRadiusPix*5), fHFRadiusPix);
mfGLF = imfilter(mfG, hGLF, 'symmetric');
mfGHF = imfilter(mfG-mfGLF, hGHF, 'symmetric');

% - Construct a cell selection image
mfCells = mfGHF;


%% -- Find cell locations by finding peaks of the cell selection image

% - Mask off the border of the image
mbAcceptMask = true(vnStackSize(1:2));
mbAcceptMask(1:ceil(fCellRadiusPix), :) = false;
mbAcceptMask(:, 1:ceil(fCellRadiusPix)) = false;
mbAcceptMask(end-ceil(fCellRadiusPix):end, :) = false;
mbAcceptMask(:, end-ceil(fCellRadiusPix):end) = false;

[mnCellLocs] = FindPeaks(mfCells, fCellRadiusPix, (mfCells > 0) & fsData.GetAlignedMask & mbAcceptMask);
nNumCells = size(mnCellLocs, 1);

% - Make some meshgrids to compute distance
vnR = 1:size(fsData, 1);
vnC = 1:size(fsData, 2);
[mnR, mnC] = ndgrid(vnR, vnC);

sCellRegions = [];
mbAlignMask = fsData.GetAlignedMask';

for (nCell = 1:nNumCells)
   mnD = sqrt((mnR - mnCellLocs(nCell, 2)).^2 + (mnC - mnCellLocs(nCell, 1)).^2);
   mbThisCellMask = (mnD <= fCellRadiusPix) & mbAlignMask;
   
   sThisRegion = bwconncomp(mbThisCellMask);
   if (isempty(sCellRegions))
      sCellRegions = sThisRegion;
   else
      sCellRegions.PixelIdxList{end+1} = find(mbThisCellMask);
      sCellRegions.NumObjects = sCellRegions.NumObjects + 1;
   end
end

% - Restore normalisation
fsData.BlankNormalisation(strNormalisation);

% - Restore warnings
warning(wOld);

return;

%% -- 

figure, imagesc(mfG'), axis equal tight, title('Raw green channel'), colorbar;
figure, imagesc(mfGLF'), axis equal tight, title('Low-pass green channel'), colorbar;
figure, imagesc(mfGHF'), axis equal tight, title('Band-pass green channel'), colorbar;



% --- END of FindCells_GRChannels.m ---
