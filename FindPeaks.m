function [mnPeakLocations] = FindPeaks(mfImage, nDiskRadius, mbAcceptMask, bWrap)

% FindPeaks - FUNCTION Find the peaks of vector or matrix, using non-maximum suppression
%
% Usage: [mnPeakLocations] = FindPeaks(mfImage, nDiskRadius <, mbAcceptMask, bWrap>)
%
% 'mfImage' is a 1D or 2D matrix.  The peaks of this matrix (over a length
% scale determined by 'nDiskRadius') will be determined.
%
% 'mbAcceptMask' is a logical matrix, the same size as 'mfImage'.  A value
% of 'true' for a pixel in 'mbAcceptMask' implies that maxima found at that
% point will be returned.
%
% 'nDiskRadius' is used to build a disk structuring element used in the
% non-maxima suppression process.  In general, maxima closer than this
% distance will be merged.
%
% 'bWrap', if true, will treat 'mfImage' as a periodic domain.  Maxima near
% the edge of 'mfImage' will be discovered in this way.  By default,
% 'bWrap' is 'false'.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: September 2008

% -- Constants

DEF_bWrap = false;


% -- Check arguments

if (~exist('mbAcceptMask', 'var') || isempty(mbAcceptMask))
   mbAcceptMask = true(size(mfImage));
end

if (~exist('bWrap', 'var') || isempty(bWrap))
   bWrap = DEF_bWrap;
end

% - How many dimensions do we have?
vnMapSize = size(mfImage);
nNumDims = sum(vnMapSize > 1);

if (nNumDims > 2)
   disp('*** FindPeaks: Only 1D and 2D input matrices are supported');
   return;
end

% - Reshape image into vector, if 1D
if (nNumDims == 1)
   mfImage = mfImage(:);
   mbAcceptMask = mbAcceptMask(:);
   vnMapSize = size(mfImage);
end

% - Wrap array, if necessary
if (bWrap)
   % - Pad by twice disk radius
   vnPadSize = nDiskRadius * ones(nNumDims, 1);
   
   % - Pad input matrix with periodic borders
   mfImagePad = padarray(mfImage, vnPadSize, 'circular');
   
   % - Pad accept mask with zeros
   mbCentralMask = padarray(mbAcceptMask, vnPadSize, 0);
else
   mfImagePad = mfImage;
   mbCentralMask = mbAcceptMask;
end

% - Make structuring element
seDisk = strel('disk', round(nDiskRadius));

% - Dilate image
mfDilated = imdilate(mfImagePad, seDisk);

% - Take central portion
mfCentral = zeros(vnMapSize);
mfCentral(mbAcceptMask) = mfDilated(mbCentralMask);

% - Locate peaks
[vYMaxima, vXMaxima] = find((abs(mfImage - mfCentral) < eps) & mbAcceptMask);

if (nNumDims == 1)
   mnPeakLocations = vYMaxima;
else
   mnPeakLocations = [vXMaxima vYMaxima];
end

% --- END of FindPeaks.m ---
