function [TIF, HEADER] = tiffread29_header(filename)

%% set defaults values :
TIF.BOS              = 'ieee-le';          %byte order string
TIF.strFilename = filename;
consolidateStrips = true;

if  isempty(findstr(filename,'.'))
   filename = [filename,'.tif'];
end

TIF.file = fopen(filename,'r','l');
if TIF.file == -1
   stkname = strrep(filename, '.tif', '.stk');
   TIF.file = fopen(stkname,'r','l');
   if TIF.file == -1
      error('tiffread2:filenotfound', ['File "',filename,'" not found.']);
   else
      filename = stkname;
   end
end
[s, m] = fileattrib(filename);

% obtain the full file path:
filename = m.Name;

%% read header
% read byte order: II = little endian, MM = big endian
byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
   TIF.BOS = 'ieee-le';                                % Intel little-endian format
elseif ( strcmp(byte_order','MM') )
   TIF.BOS = 'ieee-be';
else
   error('tiffread2:notATIF', 'This is not a TIFF file (no MM or II).');
end

%% ---- read in a number which identifies file as TIFF format
tiff_id = fread(TIF.file,1,'uint16', TIF.BOS);
if (tiff_id ~= 42)
   error('tiffread2:notATIF', 'This is not a TIFF file (missing 42).');
end

%% ---- read the byte offset for the first image file directory (IFD)
nCurrImg = 1;

HEADER(nCurrImg).img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);

nImgIndex = 1;

while  HEADER(nCurrImg).img_pos ~= 0
   
   % - Defaults
   HEADER(nCurrImg).SampleFormat     = 1;
   HEADER(nCurrImg).SamplesPerPixel  = 1;
   
   % move in the file to the first IFD
   status = fseek(TIF.file, HEADER(nCurrImg).img_pos, -1);
   if status == -1
      error('tiffread2:readerr', 'invalid file offset (error on fseek)');
   end
   
   %disp(strcat('reading img at pos :',num2str(TIF.img_pos)));
   
   %read in the number of IFD entries
   num_entries = fread(TIF.file,1,'uint16', TIF.BOS);
   %disp(strcat('num_entries =', num2str(num_entries)));
   
   %read and process each IFD entry
   for i = 1:num_entries
      
      % save the current position in the file
      file_pos  = ftell(TIF.file);
      
      % read entry tag
      entry_tag = fread(TIF.file, 1, 'uint16', TIF.BOS);
      % read entry
      entry = readIFDentry(TIF, entry_tag);
      
      switch entry_tag
         case 254
            HEADER(nCurrImg).NewSubfiletype = entry.val;
         case 256         % image width - number of column
            HEADER(nCurrImg).width          = entry.val;
         case 257         % image height - number of row
            HEADER(nCurrImg).height         = entry.val;
            HEADER(nCurrImg).ImageLength    = entry.val;
         case 258         % BitsPerSample per sample
            HEADER(nCurrImg).BitsPerSample  = entry.val;
            HEADER(nCurrImg).BytesPerSample = HEADER(nCurrImg).BitsPerSample / 8;
            HEADER(nCurrImg).bits           = HEADER(nCurrImg).BitsPerSample(1);
            %fprintf('BitsPerSample %i %i %i\n', entry.val);
         case 259         % compression
            if ( entry.val ~= 1 )
               error('tiffread2:unsupported', ['Compression format ', num2str(entry.val),' not supported.']);
            end
         case 262         % photometric interpretation
            HEADER(nCurrImg).PhotometricInterpretation = entry.val;
            if ( HEADER(nCurrImg).PhotometricInterpretation == 3 )
               warning('tiffread2:LookUp', 'Ignoring TIFF look-up table');
            end
         case 269
            HEADER(nCurrImg).document_name  = entry.val;
         case 270         % comments:
            HEADER(nCurrImg).info           = entry.val;
         case 271
            HEADER(nCurrImg).make           = entry.val;
         case 273         % strip offset
            HEADER(nCurrImg).StripOffsets   = entry.val;
            HEADER(nCurrImg).StripNumber    = entry.cnt;
            %fprintf('StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
         case 277         % sample_per pixel
            HEADER(nCurrImg).SamplesPerPixel  = entry.val;
            %fprintf('Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
         case 278         % rows per strip
            HEADER(nCurrImg).RowsPerStrip   = entry.val;
         case 279         % strip byte counts - number of bytes in each strip after any compressio
            HEADER(nCurrImg).StripByteCounts= entry.val;
         case 282         % X resolution
            HEADER(nCurrImg).x_resolution   = entry.val;
         case 283         % Y resolution
            HEADER(nCurrImg).y_resolution   = entry.val;
         case 284         %planar configuration describe the order of RGB
            HEADER(nCurrImg).PlanarConfiguration = entry.val;
         case 296         % resolution unit
            HEADER(nCurrImg).resolution_unit= entry.val;
         case 305         % software
            HEADER(nCurrImg).software       = entry.val;
         case 306         % datetime
            HEADER(nCurrImg).datetime       = entry.val;
         case 315
            HEADER(nCurrImg).artist         = entry.val;
         case 317        %predictor for compression
            if (entry.val ~= 1); error('tiffread2:unsupported', 'unsuported predictor value'); end
         case 320         % color map
            HEADER(nCurrImg).cmap           = entry.val;
            HEADER(nCurrImg).colors         = entry.cnt/3;
         case 339
            HEADER(nCurrImg).SampleFormat   = entry.val;
         case 33550       % GeoTIFF ModelPixelScaleTag
            HEADER(nCurrImg).ModelPixelScaleTag    = entry.val;
         case 33628       %metamorph specific data
            HEADER(nCurrImg).MM_private1    = entry.val;
         case 33629       %this tag identify the image as a Metamorph stack!
            TIF.MM_stack       = entry.val;
            TIF.MM_stackCnt    = entry.cnt;
         case 33630       %metamorph stack data: wavelength
            HEADER(nCurrImg).MM_wavelength  = entry.val;
         case 33631       %metamorph stack data: gain/background?
            HEADER(nCurrImg).MM_private2    = entry.val;
         case 33922       % GeoTIFF ModelTiePointTag
            IMG.ModelTiePointTag    = entry.val;
         case 34412       % Zeiss LSM data
            HEADER(nCurrImg).LSM_info           = entry.val;
         case 34735       % GeoTIFF GeoKeyDirectory
            HEADER(nCurrImg).GeoKeyDirTag       = entry.val;
         case 34737       % GeoTIFF GeoASCIIParameters
            HEADER(nCurrImg).GeoASCII       = entry.val;
         case 42113       % GeoTIFF GDAL_NODATA
            HEADER(nCurrImg).GDAL_NODATA    = entry.val;
         otherwise
            warning('tiffread2:ignored_tag', 'Ignored TIFF entry with tag %i (cnt %i)\n', entry_tag, entry.cnt);
      end
      
      % calculate bounding box  if you've got the stuff
      if isfield(HEADER(nCurrImg), 'ModelPixelScaleTag') && isfield(HEADER(nCurrImg), 'ModelTiePointTag') && isfield(HEADER(nCurrImg), 'height')&& isfield(HEADER(nCurrImg), 'width'),
         HEADER(nCurrImg).North=HEADER(nCurrImg).ModelTiePointTag(5)-HEADER(nCurrImg).ModelPixelScaleTag(2)*HEADER(nCurrImg).ModelTiePointTag(2);
         HEADER(nCurrImg).South=HEADER(nCurrImg).North-HEADER(nCurrImg).height*HEADER(nCurrImg).ModelPixelScaleTag(2);
         HEADER(nCurrImg).West=HEADER(nCurrImg).ModelTiePointTag(4)+HEADER(nCurrImg).ModelPixelScaleTag(1)*HEADER(nCurrImg).ModelTiePointTag(1);
         HEADER(nCurrImg).East=HEADER(nCurrImg).West+HEADER(nCurrImg).width*HEADER(nCurrImg).ModelPixelScaleTag(1);
      end
      
      % move to next IFD entry in the file
      status = fseek(TIF.file, file_pos+12, -1);
      if status == -1
         error('tiffread2:readerr', 'invalid file offset (error on fseek)');
      end
   end
   
   %Planar configuration is not fully supported
   %Per tiff spec 6.0 PlanarConfiguration irrelevent if SamplesPerPixel==1
   %Contributed by Stephen Lang
   if (HEADER(nCurrImg).SamplesPerPixel ~= 1) && ( ~isfield(HEADER(nCurrImg), 'PlanarConfiguration') || HEADER(nCurrImg).PlanarConfiguration == 1 )
      error('tiffread2:unsupported', 'PlanarConfiguration = 1 is not supported');
   end
   
   %total number of bytes per image:
   PlaneBytesCnt = HEADER(nCurrImg).width * HEADER(nCurrImg).height * HEADER(nCurrImg).BytesPerSample;
   
   %% try to consolidate the TIFF strips if possible
   
   if consolidateStrips
      %Try to consolidate the strips into a single one to speed-up reading:
      BytesCnt = HEADER(nCurrImg).StripByteCounts(1);
      
      if BytesCnt < PlaneBytesCnt
         
         ConsolidateCnt = 1;
         %Count how many Strip are needed to produce a plane
         while HEADER(nCurrImg).StripOffsets(1) + BytesCnt == HEADER(nCurrImg).StripOffsets(ConsolidateCnt+1)
            ConsolidateCnt = ConsolidateCnt + 1;
            BytesCnt = BytesCnt + HEADER(nCurrImg).StripByteCounts(ConsolidateCnt);
            if ( BytesCnt >= PlaneBytesCnt ); break; end
         end
         
         %Consolidate the Strips
         if ( BytesCnt <= PlaneBytesCnt(1) ) && ( ConsolidateCnt > 1 )
            %fprintf('Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
            HEADER(nCurrImg).StripByteCounts = [BytesCnt; HEADER(nCurrImg).StripByteCounts(ConsolidateCnt+1:HEADER(nCurrImg).StripNumber ) ];
            HEADER(nCurrImg).StripOffsets = HEADER(nCurrImg).StripOffsets( [1 , ConsolidateCnt+1:HEADER(nCurrImg).StripNumber] );
            HEADER(nCurrImg).StripNumber  = 1 + HEADER(nCurrImg).StripNumber - ConsolidateCnt;
         end
      end
   end
   
   
   %% Record the next IFD address:
   HEADER(nCurrImg+1).img_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);
   %if (TIF.img_pos) disp(['next ifd at', num2str(TIF.img_pos)]); end
   
   % - Record the start of the image data for this image
   HEADER(nCurrImg).nImageOffset = ftell(TIF.file);
   
   % - Skip to next image
   nCurrImg = nCurrImg + 1;
end


% clean-up
HEADER = HEADER(1:end-1);

% fclose(TIF.file);

end



%% ==================sub-functions that reads an IFD entry:===================


function [nbBytes, matlabType] = convertType(tiffType)
switch (tiffType)
   case 1
      nbBytes=1;
      matlabType='uint8';
   case 2
      nbBytes=1;
      matlabType='uchar';
   case 3
      nbBytes=2;
      matlabType='uint16';
   case 4
      nbBytes=4;
      matlabType='uint32';
   case 5
      nbBytes=8;
      matlabType='uint32';
   case 7
      nbBytes=1;
      matlabType='uchar';
   case 11
      nbBytes=4;
      matlabType='float32';
   case 12
      nbBytes=8;
      matlabType='float64';
   otherwise
      error('tiffread2:unsupported', 'tiff type %i not supported', tiffType)
end
end

%% ==================sub-functions that reads an IFD entry:===================

function  entry = readIFDentry(TIF, entry_tag)

entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.BOS);
entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.BOS);
%disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);

[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
   %next field contains an offset:
   offset = fread(TIF.file, 1, 'uint32', TIF.BOS);
   %disp(strcat('offset = ', num2str(offset)));
   status = fseek(TIF.file, offset, -1);
   if status == -1
      error('tiffread2:corrupted', 'invalid file offset (error on fseek)');
   end
   
end


if entry_tag == 33629   % metamorph 'rationals'
   entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.BOS);
elseif entry_tag == 34412  %TIF_CZ_LSMINFO
   entry.val = readLSMinfo(TIF);
else
   if entry.tiffType == 5
      entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);
   else
      entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.BOS);
   end
end

if ( entry.tiffType == 2 );
   entry.val = char(entry.val');
end

end




%% ==============partial-parse of LSM info:

function R = readLSMinfo(TIF)

% Read part of the LSM info table version 2
% this provides only very partial information, since the offset indicate that
% additional data is stored in the file

R.MagicNumber            = sprintf('0x%09X',fread(TIF.file, 1, 'uint32', TIF.BOS));
StructureSize          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionX             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionY             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionZ             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionChannels      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionTime          = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.IntensityDataType      = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ThumbnailX             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ThumbnailY             = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.VoxelSizeX             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeY             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeZ             = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginX                = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginY                = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginZ                = fread(TIF.file, 1, 'float64', TIF.BOS);
R.ScanType               = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.SpectralScan           = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.DataType               = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetVectorOverlay    = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetInputLut         = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetOutputLut        = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetChannelColors    = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.TimeInterval           = fread(TIF.file, 1, 'float64', TIF.BOS);
OffsetChannelDataTypes = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetScanInformation  = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetKsData           = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetTimeStamps       = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetEventList        = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetRoi              = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetBleachRoi        = fread(TIF.file, 1, 'uint32', TIF.BOS);
OffsetNextRecording    = fread(TIF.file, 1, 'uint32', TIF.BOS);

% There are more information stored in this table, which is not read here


%read real acquisition times:
if ( OffsetTimeStamps > 0 )
   
   status = fseek(TIF.file, OffsetTimeStamps, -1);
   if status == -1
      error('tiffread2:readerror', 'error on fseek');
   end
   
   StructureSize          = fread(TIF.file, 1, 'int32', TIF.BOS);
   NumberTimeStamps       = fread(TIF.file, 1, 'int32', TIF.BOS);
   for i=1:NumberTimeStamps
      R.TimeStamp(i)       = fread(TIF.file, 1, 'float64', TIF.BOS);
   end
   
   %calculate elapsed time from first acquisition:
   R.TimeOffset = R.TimeStamp - R.TimeStamp(1);
   
end


end
