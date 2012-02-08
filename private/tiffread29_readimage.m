function [data, stack] = tiffread29_readimage(TIF, HEADER, vnFrames)

% - Preallocate data block
classname = get_data_class(HEADER(vnFrames(1)), 1);
data = zeros(HEADER(vnFrames(1)).height, HEADER(vnFrames(1)).width, numel(vnFrames), HEADER(vnFrames(1)).SamplesPerPixel, classname);
stack = [];

for (nFrame = numel(vnFrames):-1:1) % Go backwards to pre-allocate
   % - Skip to image offset
   fseek(TIF.file, HEADER(nFrame).nImageOffset, 'bof');
   
   if isfield( TIF, 'MM_stack' )
      % This part reads a metamorph TIF
      
      % Only take existing frames
      vnFrames = vnFrames(vnFrames <= TIF.MM_stackCnt);
      
      %this loop reads metamorph stacks:
      for ii = vnFrames
         
         StripCnt = 1;
         offset = PlaneBytesCnt * (ii-1);
         
         %read the image channels
         for c = 1:HEADER(vnFrames(nFrame)).SamplesPerPixel
            data(:, :, nFrame, c) = read_plane(TIF, HEADER(vnFrames(nFrame)), offset, c, StripCnt);
         end
         
         [ IMG.MM_stack, IMG.MM_wavelength, IMG.MM_private2 ] = splitMetamorph(ii, TIF);
         
         stack(nFrame) = IMG; %#ok<AGROW>
         
      end
      
   else
      
      %this part reads a normal TIFF stack:
      
      StripCnt = 1;
      %read the image channels
      for c = 1:HEADER(nFrame).SamplesPerPixel
         data(:, :, nFrame, c) = read_plane(TIF, HEADER(vnFrames(nFrame)), 0, c, StripCnt);
      end
   end
end

end

%% ===========================================================================

function plane = read_plane(TIF, HEADER_frame, offset, plane_nb, StripCnt)

%return an empty array if the sample format has zero bits
if ( HEADER_frame.BitsPerSample(plane_nb) == 0 )
   plane=[];
   return;
end

%fprintf('reading plane %i size %i %i\n', plane_nb, width, height);

classname = get_data_class(HEADER_frame, plane_nb);

% Preallocate a matrix to hold the sample data:
try
   plane = zeros(HEADER_frame.height, HEADER_frame.width, classname);
catch
   %compatibility with older matlab versions:
   eval(['plane = ', classname, '(zeros(HEADER_frame.height, HEADER_frame.width));']);
end

% Read the strips and concatenate them:
line = 1;
while ( StripCnt <= HEADER_frame.StripNumber )
   
   strip = read_strip(TIF, HEADER_frame, offset, plane_nb, StripCnt, classname);
   StripCnt = StripCnt + 1;
   
   % copy the strip onto the data
   plane(line:(line+size(strip,1)-1), :) = strip;
   
   line = line + size(strip,2);
   if ( line > HEADER_frame.height )
      break;
   end
   
end

% Extract valid part of data if needed
if ~all(size(plane) == [HEADER_frame.height HEADER_frame.width]),
   plane = plane(1:HEADER_frame.height, 1:HEADER_frame.width);
   warning('tiffread2:Crop','Cropping data: found more bytes than needed');
end

end


%% ================== sub-functions to read a strip ===================

function strip = read_strip(TIF, HEADER_frame, offset, plane_nb, stripCnt, classname)

%fprintf('reading strip at position %i\n',TIF.StripOffsets(stripCnt) + offset);
StripLength = HEADER_frame.StripByteCounts(stripCnt) ./ HEADER_frame.BytesPerSample(plane_nb);

%fprintf( 'reading strip %i\n', stripCnt);
status = fseek(TIF.file, HEADER_frame.StripOffsets(stripCnt) + offset, 'bof');
if status == -1
   error('tiffread2:readerr', 'invalid file offset (error on fseek)');
end

bytes = fread( TIF.file, StripLength, classname, TIF.BOS );

if any( length(bytes) ~= StripLength )
   error('tiffread2:corrupted', 'End of file reached unexpectedly.');
end

strip = reshape(bytes, HEADER_frame.width, StripLength / HEADER_frame.width)';

end

%% =============distribute the metamorph infos to each frame:
function [MMstack, MMwavelength, MMprivate2] = splitMetamorph(imgCnt)

global TIF;

MMstack = [];
MMwavelength = [];
MMprivate2 = [];

if TIF.MM_stackCnt == 1
   return;
end

left  = imgCnt - 1;

if isfield( TIF, 'MM_stack' )
   S = length(TIF.MM_stack) / TIF.MM_stackCnt;
   MMstack = TIF.MM_stack(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_wavelength' )
   S = length(TIF.MM_wavelength) / TIF.MM_stackCnt;
   MMwavelength = TIF.MM_wavelength(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_private2' )
   S = length(TIF.MM_private2) / TIF.MM_stackCnt;
   MMprivate2 = TIF.MM_private2(S*left+1:S*left+S);
end

end


%% === determine the data class for this frame

function strClass = get_data_class(HEADER_frame, plane_nb)

%determine the type needed to store the pixel values:
switch( HEADER_frame.SampleFormat )
   case 1
      strClass = sprintf('uint%i', HEADER_frame.BitsPerSample(plane_nb));
   case 2
      strClass = sprintf('int%i', HEADER_frame.BitsPerSample(plane_nb));
   case 3
      if ( HEADER_frame.BitsPerSample(plane_nb) == 32 )
         strClass = 'single';
      else
         strClass = 'double';
      end
   otherwise
      error('tiffread2:unsupported', 'unsuported TIFF sample format %i', HEADER_frame.SampleFormat);
end

end
