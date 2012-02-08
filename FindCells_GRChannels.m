function [sCellRegions] = FindCells_GRChannels(fsData, fLFRadiusUM, fHFRadiusUM, fCellRadiusUM)

% FindCells_GRChannels - FUNCTION Try to locate neuron somata using red-channel subtraction and filtering
%
% Usage: [sCellRegions] = FindCells_GRChannels(fsData <, fLFRadiusUM, fHFRadiusUM, fCellRadiusUM>)
%
% 'fsData' is a FocusStack object.  The red channel (assumed to label astrocytes
% and not label neurons) will be subtracted from the green channel (assumed to
% label everything).  The average images (taken from stimulus response frames
% only, if the stimulation times are available in the stack) will be used.
% Band-pass filtering will be used to locate neurons.
%
% 'sCellRegions' will be a region structure, as returned by matlab function
% bwconncomp, defining the circular regions corresponding to each cell.
%
% The optional parameters 'fLFRadiusUM', 'fHFRadiusUM' and 'fCellRadiusUM'
% control the filtering process.  By default, they adopt reasonable values for
% sulfarhodamine plus injected OGB labelling in mouse visual cortex.  All values
% are in micrometres, and rely on the spatial calibration of the stack.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: January 2011

% -- Default parameters

DEF_fCellRadiusUM = 2.5;
DEF_fHFRadiusUM = 1.6;
DEF_fLFRadiusUM = 10;


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

% - Get average green and red channel responses
mfG = double(fsData.SummedAlignedFrames(:, :, vbUseFrame, 1)) ./ vnStackSize(3);
mfR = double(fsData.SummedAlignedFrames(:, :, :, 2)) ./ vnStackSize(3);

% - Construct filters, and filter average frames
hGLF = fspecial('gaussian', ceil(fLFRadiusPix*5), fLFRadiusPix);
hGHF = fspecial('gaussian', ceil(fLFRadiusPix*5), fHFRadiusPix);

mfGLF = imfilter(mfG, hGLF, 'symmetric');
mfRLF = imfilter(mfR, hGLF, 'symmetric');

mfGHF = imfilter(mfG-mfGLF, hGHF, 'symmetric');
mfRHF = imfilter(mfR-mfRLF, hGHF, 'symmetric');

% mfGHF = mfGHF - mean(mfGHF(:));
% mfGHF(mfGHF < 0) = 0;
% mfGHF = mfGHF ./ max(mfGHF(:));
% 
% mfRHF = -mfRHF + mean(mfRHF(:));
% % mfRHF(mfRHF < 0) = 0;
% mfRHF = mfRHF ./ max(mfRHF(:));

% - Construct a cell selection image
mfCells = (mfGHF) - (mfRHF);% + (mfGHF .* mfRHF);


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
mbAlignMask = fsData.GetAlignedMask;

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

% return;

%% -- 

figure, imagesc(mfG'), axis equal tight, title('Raw green channel'), colorbar;
figure, imagesc(mfGLF'), axis equal tight, title('Low-pass green channel'), colorbar;
figure, imagesc(mfGHF'), axis equal tight, title('Band-pass green channel'), colorbar;

figure, imagesc(mfR'), axis equal tight, title('Raw red channel'), colorbar;
figure, imagesc(mfRLF'), axis equal tight, title('Low-pass red channel'), colorbar;
figure, imagesc(mfRHF'), axis equal tight, title('Band-pass red channel'), colorbar;


im = (mfR-mfRLF)';
% im(im < median(im(:))/2) = median(im(:));
im = im - min(im(:));
im = im ./ max(im(:));
mfGC = (mfG-mfGLF)';
% mfGC(mfGC < median(mfGC(:))) = median(mfGC(:));
mfGC = mfGC - min(mfGC(:));
mfGC = mfGC ./ max(mfGC(:));
im(:, :, 2) = mfGC;
im(:, :, 3) = 0;

figure, image(im);
hold on;
axis equal tight off;
plot(mnCellLocs(:, 2), mnCellLocs(:, 1), 'w.');

figure, imagesc(mfCells'), axis equal tight, title('Criteria image'), colorbar;
hold on;
plot(mnCellLocs(:, 2), mnCellLocs(:, 1), 'w.');


% --- END of FindCells_GRChannels.m ---
