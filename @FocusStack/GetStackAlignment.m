function [mfFrameOffsets] = GetStackAlignment(tfStack, vnChannel, bProgressive, nUpsampling, mfReferenceImage)

% GetStackAlignment - METHOD Compute registration for an image stack (possible sub-pixel resolution)
%
% Usage: [mfFrameOffsets] = GetStackAlignment(tfStack <, vnChannel, bProgressive, nUpsampling, mfReferenceImage>)
%
% 'tfStack' is either a tensor image stack [X Y F] or [X Y F C], where F is a
% frame index and C is an optional channel index.
%
% The optional argument 'vnChannel' specifies which stack channel to use for
% registration.  If you want to use several channels summed, specify a vector of
% channel indices.
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
% than the initial stack frame.
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

if (exist('mfReferenceImage', 'var'))
   if (bProgressive)
      error('FocusStack:AlignInvalidArguments', '*** FocusStack/GetStackAlignment: Cannot provide a reference image for progressive alignment.');
   end
   
   if (~isequal(size(mfReferenceImage), vnStackSize(1:2)))
      error('FocusStack:InvalidReference', '*** FocusStack/GetStackAlignment: Reference image is the wrong size.');
   end
end


% -- Turn off data normalisation

if (isa(tfStack, 'FocusStack'))
   bSubtractBlack = tfStack.bSubtractBlack;
   strNormalisation = tfStack.BlankNormalisation('none');
end


% -- Determine registration with dftregistration, for each frame

nNumFrames = size(tfStack, 3);
mfFrameOffsets = zeros(nNumFrames, 2);

% - Compute first frame FFT, if doing first frame only registration
if (~bProgressive)
   if (~exist('mfReferenceImage', 'var'))
      mfFFTInitialFrame = fft2(sum(tfStack.ExtractFrames({':', ':', 1, vnChannel}), 4));
   else
      mfFFTInitialFrame = fft2(mfReferenceImage);
   end
end

% - Show some progress
fprintf(1, 'Measuring stack alignment: %3d%%', 0);

for (nFrame = 1:nNumFrames)
   % - Compute the FFT for this frame
   mfFFTThisFrame = fft2(sum(tfStack.ExtractFrames({':', ':', nFrame, vnChannel}), 4));

   if (bProgressive)
      % - Register the last frame against this one
      mfFFTRegFrame = fft2(sum(tfStack.ExtractFrames({':', ':', nFrame-1, vnChannel}), 4)); %#ok<PFBNS>
   else
      % - Register this frame against the first frame
      mfFFTRegFrame = mfFFTInitialFrame;
   end
   
   % - Determine the registration
   [vOutput] = dftregistration(mfFFTRegFrame, mfFFTThisFrame, nUpsampling);
   mfFrameOffsets(nFrame, :) = vOutput([4 3]);
   
   % - Show some progress
   fprintf(1, '\b\b\b\b%3d%%', round(nFrame / nNumFrames * 100));
end

fprintf(1, '\n');

% - Convert to global offsets, if required 
if (bProgressive)
   mfFrameOffsets = cumsum(mfFrameOffsets, 1);
end


% -- Restore data normalisation

if (isa(tfStack, 'FocusStack'))
    tfStack.BlankNormalisation(strNormalisation);
end

% --- END of GetStackAlignemnt METHOD ---
