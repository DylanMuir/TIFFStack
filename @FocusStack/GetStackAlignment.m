function [mfFrameOffsets] = GetStackAlignment(tfStack, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPUM)

% GetStackAlignment - METHOD Compute registration for an image stack (possible sub-pixel resolution)
%
% Usage: [mfFrameOffsets] = GetStackAlignment(tfStack <, vnChannel, bProgressive, nUpsampling, mfReferenceImage, nWindowLength, vfSpatFreqCutoffCPP>)
%        [mfFrameOffsets] = GetStackAlignment(tfStack <, vnChannel, bProgressive, nUpsampling, nReferenceFrame, nWindowLength, vfSpatFreqCutoffCPP>)
%        [mfFrameOffsets] = GetStackAlignment(tfStack <, cvnChannelComb, ...)
%
% 'tfStack' is either a tensor image stack [X Y F] or [X Y F C], where F is a
% frame index and C is an optional channel index.
%
% The optional argument 'vnChannel' specifies which stack channel to use for
% registration.  If you want to use several channels summed, specify a vector of
% channel indices.
%
% 'cvnChannelComb' can be provided instead; this is a cell array, where the
% first cell contains a function handle that will be applied to a tensor [X Y F
% C], and must produce an output [X Y F] by combining channels in some way.  The
% second cell contains a vector of channel indices to use.  For example:
%    {@(t)nansum(t, 4) [1 2]}
% would extract channels 1 and 2 then sum them, and perform the alignment on
% the resulting image.
%
% The optional argument 'bProgressive' determines whether shifts should be
% calculated between successive frames ('bProgressive' = true), or whether all
% shifts should be calculated with the first frame ('bProgressive' = false,
% default).
%
% The optional argument 'nUpsampling' specifies that the registration should be
% computed to 1/'nUpsampling' pixels.  Default is 1, meaning that registration
% occurs to single pixel resolution.
%
% 'mfReferenceImage' is an optional reference image used for alignment, rather
% than the initial stack frame.  Optionally, a scalar frame index can be
% supplied instead.  In this case, the indicated frame will be extracted from
% the stack to be used as a reference for alignment.
%
% 'nWindowLength' is an optional parameter that specifies the number of frames
% to average together, in a moving window, to determine the alignment for the
% current frame.  Default is 1, meaning that only a single frame is used.
%
% 'vfSpatFreqCutoffCPP' is an optional parameter that defines the cutoff spatial
% frequencies to include in computing the mis-alignment.  The vector is
% [fMinFreqCPP fMaxFreqCPP], both in cycles per pixel.  Spatial frequencies
% between these limits are included by a band-pass filter.
%
% 'mfFrameOffsets' gives the offsets [x y] to shift each frame, such that the
% whole stack is in alignment.
%
% Uses 'dftregistration'.
%
% Citation: Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 8th November, 2010

% -- Defaults

DEF_vnChannel = 1;
DEF_bProgressive = false;
DEF_nUpsampling = 1;
DEF_vfSpatFreqCutoffCPUM = [0 inf];


% -- Check arguments

vnStackSize = size(tfStack);

if (nargin < 1)
   disp('*** GetStackAlignment: Incorrect usage');
   help GetStackAlignment;
   return;
end

if (~exist('vnChannel', 'var') || isempty(vnChannel))
   vnChannel = DEF_vnChannel;
end

if (~exist('bProgressive', 'var') || isempty(bProgressive))
   bProgressive = DEF_bProgressive;
end

if (~exist('nUpsampling', 'var') || isempty(nUpsampling))
   nUpsampling = DEF_nUpsampling;
end

if (iscell(vnChannel))
   % - An arbitrary channel combination function was supplied
   fhChannelFunc = vnChannel{1};
   vnChannel = vnChannel{2};
   
elseif (isscalar(vnChannel))
   % - Single channel, so just return the frame
   fhChannelFunc = @(t)t;

else
   % - Several channels, so return the sum
   fhChannelFunc = @(t)nansum(t, 4);
end


% -- Turn off data normalisation

if (isa(tfStack, 'FocusStack'))
   bSubtractBlack = tfStack.bSubtractBlack;
   strNormalisation = tfStack.BlankNormalisation('none');
end


% -- Make FFT frequency mesh

fSpatSamplingPPUM = 1;
vnStimSize = size(tfStack);
vnSFiltSize = 2.^ [nextpow2(vnStimSize(1)) nextpow2(vnStimSize(2)) nextpow2(vnStimSize(3))];
vfxSFreq = ifftshift(fSpatSamplingPPUM .* ((-vnSFiltSize(1)/2) : (vnSFiltSize(1)/2-1)) ./ vnSFiltSize(1));
vfySFreq = ifftshift(fSpatSamplingPPUM .* ((-vnSFiltSize(2)/2) : (vnSFiltSize(2)/2-1)) ./ vnSFiltSize(2));
[mfYSFreq, mfXSFreq] = meshgrid(vfySFreq, vfxSFreq);

% - Generate a spatial band-pass filter
mfSpatFilter = false(vnSFiltSize(1:2));
mfSpatFilter(sqrt((mfXSFreq.^2) + (mfYSFreq.^2)) >= vfSpatFreqCutoffCPUM(1)) = 1;
mfSpatFilter(sqrt((mfXSFreq.^2) + (mfYSFreq.^2)) > vfSpatFreqCutoffCPUM(2)) = 0;


% -- Determine registration with dftregistration, for each frame

nNumFrames = size(tfStack, 3);
mfFrameOffsets = zeros(nNumFrames, 2);

% - Determine frame window for computing alignment
vnWindow = ceil(-(nWindowLength-1)/2):floor(nWindowLength/2);
nFirstFrame = 1-vnWindow(1);
nLastFrame = size(tfStack, 3) - vnWindow(end);

% - Make initial window
tfWindow = fhChannelFunc(tfStack.ExtractFrames({':', ':', vnWindow + nFirstFrame, vnChannel}));
nWindowIndex = 1;

if (exist('mfReferenceImage', 'var') && ~isempty(mfReferenceImage))
   if (isscalar(mfReferenceImage) && (mfReferenceImage > 0) && (mfReferenceImage < size(tfStack, 3)))
      % - Assume it's a frame index, so extract the appropriate frame(s)
      vnRefWindow = vnWindow + mfReferenceImage;
      vnRefWindow = vnRefWindow(vnRefWindow >= 1);
      vnRefWindow = vnRefWindow(vnRefWindow <= size(tfStack, 3));
      mfReferenceImage = nansum(fhChannelFunc(tfStack.ExtractFrames({':', ':', vnRefWindow, vnChannel})), 3);
      
   else
      if (bProgressive)
         error('FocusStack:AlignInvalidArguments', '*** FocusStack/GetStackAlignment: Cannot provide a reference image for progressive alignment.');
      end
      
      if (~isequal(size(mfReferenceImage), vnStackSize(1:2)))
         error('FocusStack:InvalidReference', '*** FocusStack/GetStackAlignment: Reference image is the wrong size.');
      end
   end
end


% - Compute first frame FFT
if (exist('mfReferenceImage', 'var') && ~isempty(mfReferenceImage))
   mfFFTRegFrame = fft2(mfReferenceImage);
else
   mfFFTRegFrame = fft2(nansum(tfWindow, 3));
end

mfFFTRegFrame = mfFFTRegFrame .* mfSpatFilter;

% - Show some progress
fprintf(1, 'Measuring stack alignment: %3d%%', 0);

% hFig = figure;

for (nFrame = (nFirstFrame+1):nLastFrame)
   % - Build up window
   tfWindow(:, :, nWindowIndex) = fhChannelFunc(tfStack.ExtractFrames({':', ':', nFrame, vnChannel}));
   
   % - Compute the FFT for this window
   mfFFTThisFrame = fft2(nansum(tfWindow, 3));

   % - Perform a spatial band-pass filter on the frame
   mfFFTThisFrame = mfFFTThisFrame .* mfSpatFilter;
   
%    figure(hFig);
%    imagesc(ifft2(mfFFTThisFrame));
%    drawnow;
   
   % - Determine the registration
   [vOutput] = dftregistration(mfFFTRegFrame, mfFFTThisFrame, nUpsampling);
   mfFrameOffsets(nFrame, :) = vOutput([4 3]);

   % - Record registration frame, if doing progressive alignment
   if (bProgressive)
      % - Register the next frame against this one
      mfFFTRegFrame = mfFFTThisFrame;
   end

   % - Move to next window index
   nWindowIndex = mod(nWindowIndex, nWindowLength) + 1;
   
   % - Show some progress
   fprintf(1, '\b\b\b\b%3d%%', round(nFrame / nNumFrames * 100));
end

fprintf(1, '\n');

% - Convert to global offsets, if required 
if (bProgressive)
   mfFrameOffsets = cumsum(mfFrameOffsets, 1);
end

% - Align initial and final frames
mfFrameOffsets(1:(nFirstFrame-1), 1) = mfFrameOffsets(nFirstFrame, 1);

mfFrameOffsets(nLastFrame+1:end, 1) = mfFrameOffsets(nLastFrame, 1);
mfFrameOffsets(nLastFrame+1:end, 2) = mfFrameOffsets(nLastFrame, 2);

% -- Restore data normalisation

if (isa(tfStack, 'FocusStack'))
    tfStack.BlankNormalisation(strNormalisation);
end

% --- END of GetStackAlignemnt METHOD ---
