function [TIF, HEADER] = tiffread31_header(file_name)

%% set defaults values :

opt.ReadUnknownTags = true;
opt.ConsolidateStrips = true;
opt.SimilarImages = false;
opt.DistributeMetaData = true;

% the structure IMG is returned to the user, while TIF is not.
% so tags usefull to the user should be stored as fields in IMG, while
% those used only internally can be stored in TIF.
% the structure ANDOR has additional header information which is added to
% each plane of the image eventually

TIF = struct('ByteOrder', 'ieee-le');          %byte order string
TIF.SamplesPerPixel = 1;
TIF.PlanarConfiguration = 1;

% obtain the full file path:
[status, file_attrib] = fileattrib(file_name);

if status == 0
   error('tiffread3:filenotfound', ['File "',file_name,'" not found.']);
end

file_name = file_attrib.Name;

% open file for reading
TIF.file = fopen(file_name,'r','l');

% obtain the short file name:
[~, name, ext] = fileparts(file_name);
image_name = [name, ext];


%% ---- read byte order: II = little endian, MM = big endian

bos = fread(TIF.file, 2, '*char');
if ( strcmp(bos', 'II') )
   TIF.ByteOrder = 'ieee-le';    % Intel little-endian format
elseif ( strcmp(bos','MM') )
   TIF.ByteOrder = 'ieee-be';
else
   error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies TIFF format

tiff_id = fread(TIF.file,1,'uint16', TIF.ByteOrder);

if (tiff_id ~= 42)
   error('This is not a TIFF file (missing 42).');
end

%% ---- read the image file directories (IFDs)

ifd_pos   = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
img_indx  = 0;

while  ifd_pos ~= 0
   
   img_indx = img_indx + 1;
   
   TIF.file_name = file_name; %#ok<*AGROW>
   TIF.image_name = image_name;
   HEADER(img_indx).index = img_indx;
   
   % move in the file to the next IFD
   file_seek(ifd_pos);
   
   %read in the number of IFD entries
   num_entries = fread(TIF.file,1,'uint16', TIF.ByteOrder);
   %fprintf('num_entries = %i\n', num_entries);
   
   % store current position:
   entry_pos = ifd_pos+2;
   
   % read the next IFD address:
   file_seek(ifd_pos+12*num_entries+2);
   ifd_pos = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
   
   %fprintf('reading IFD %i at position %i\n', img_indx, ifd_pos);
   
   % read all the IFD entries
   for inx = 1:num_entries
      
      % move to next IFD entry in the file
      file_seek(entry_pos+12*(inx-1));
      
      % read entry tag
      entry_tag = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
      % read entry
      entry = readIFDentry(entry_tag);
      
      %fprintf('found entry with tag %i\n', entry_tag);
      
      
      % not all valid tiff tags have been included, but tags can easily be added to this code
      % See the official list of tags:
      % http://partners.adobe.com/asn/developer/pdfs/tn/TIFF6.pdf
      switch entry_tag
         
         case 254
            TIF.NewSubfiletype = entry.val;
         case 256         % image width = number of column
            HEADER(img_indx).width = entry.val;
         case 257         % image height = number of row
            HEADER(img_indx).height = entry.val;
            TIF.ImageLength = entry.val;
         case 258
            TIF.BitsPerSample = entry.val;
            TIF.BytesPerSample = TIF.BitsPerSample / 8;
            HEADER(img_indx).bits = TIF.BitsPerSample(1);
            %fprintf('BitsPerSample %i %i %i\n', entry.val);
         case 259         % compression
            if ( entry.val ~= 1 )
               error(['Compression format ', num2str(entry.val),' not supported.']);
            end
         case 262         % photometric interpretation
            TIF.PhotometricInterpretation = entry.val;
            if ( TIF.PhotometricInterpretation == 3 )
               %warning('tiffread:LookUp', 'Ignoring TIFF look-up table');
               fprintf(2, 'Ignoring TIFF look-up table in %s\n', image_name);
            end
         case 269
            HEADER(img_indx).document_name = entry.val;
         case 270         % general comments:
            HEADER(img_indx).info = entry.val;
         case 271
            HEADER(img_indx).make = entry.val;
         case 273         % strip offset
            HEADER(img_indx).StripOffsets = entry.val;
            HEADER(img_indx).StripNumber = entry.cnt;
            %fprintf('StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
         case 274
            % orientation is read, but the matrix is not rotated
            if ( 1 < entry.val ) && ( entry.val < 9 )
               HEADER(img_indx).orientation = entry.val;
               keys = {'TopLeft', 'TopRight', 'BottomRight', 'BottomLeft', 'LeftTop', 'RightTop', 'RightBottom', 'LeftBottom'};
               HEADER(img_indx).orientation_text = keys{entry.val};
            end
         case 277
            TIF.SamplesPerPixel  = entry.val;
            %fprintf('Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
         case 278
            HEADER(img_indx).RowsPerStrip   = entry.val;
         case 279         % strip byte counts - number of bytes in each strip after any compressio
            HEADER(img_indx).StripByteCounts= entry.val;
         case 282
            HEADER(img_indx).x_resolution   = entry.val;
         case 283
            HEADER(img_indx).y_resolution   = entry.val;
         case 284         %planar configuration describe the order of RGB
            TIF.PlanarConfiguration = entry.val;
         case 296
            HEADER(img_indx).resolution_unit= entry.val;
         case 305
            HEADER(img_indx).software       = entry.val;
         case 306
            HEADER(img_indx).datetime       = entry.val;
         case 315
            HEADER(img_indx).artist         = entry.val;
         case 317        %predictor for compression
            if (entry.val ~= 1); error('unsuported predictor value'); end
         case 320         % color map
            HEADER(img_indx).cmap           = entry.val;
            HEADER(img_indx).colors         = entry.cnt/3;
         case 339
            TIF.SampleFormat   = entry.val;
            
            % ANDOR tags
         case 4864
            if ~exist('ANDOR', 'var')
               %the existence of the variable indicates that we are
               %handling a Andor generated file
               HEADER(img_indx).ANDOR = struct([]);
            end
         case 4869       %ANDOR tag: temperature in Celsius when stabilized
            if ~(entry.val == -999)
               HEADER(img_indx).ANDOR.temperature = entry.val;
            end
         case 4876       %exposure time in seconds
            HEADER(img_indx).ANDOR.exposureTime   = entry.val;
         case 4878
            HEADER(img_indx).ANDOR.kineticCycleTime = entry.val;
         case 4879       %number of accumulations
            HEADER(img_indx).ANDOR.nAccumulations = entry.val;
         case 4881
            HEADER(img_indx).ANDOR.acquisitionCycleTime = entry.val;
         case 4882       %Readout time in seconds, 1/readoutrate
            HEADER(img_indx).ANDOR.readoutTime = entry.val;
         case 4884
            if (entry.val == 9)
               HEADER(img_indx).ANDOR.isPhotonCounting = 1;
            else
               HEADER(img_indx).ANDOR.isPhotonCounting = 0;
            end
         case 4885         %EM DAC level
            HEADER(img_indx).ANDOR.emDacLevel = entry.val;
         case 4890
            HEADER(img_indx).ANDOR.nFrames = entry.val;
         case 4896
            HEADER(img_indx).ANDOR.isFlippedHorizontally = entry.val;
         case 4897
            HEADER(img_indx).ANDOR.isFlippedVertically = entry.val;
         case 4898
            HEADER(img_indx).ANDOR.isRotatedClockwise = entry.val;
         case 4899
            HEADER(img_indx).ANDOR.isRotatedAnticlockwise = entry.val;
         case 4904
            HEADER(img_indx).ANDOR.verticalClockVoltageAmplitude = entry.val;
         case 4905
            HEADER(img_indx).ANDOR.verticalShiftSpeed = entry.val;
         case 4907
            HEADER(img_indx).ANDOR.preAmpSetting = entry.val;
         case 4908         %Camera Serial Number
            HEADER(img_indx).ANDOR.serialNumber = entry.val;
         case 4911       %Actual camera temperature when not equal to -999
            if ~(entry.val == -999)
               HEADER(img_indx).ANDOR.unstabilizedTemperature = entry.val;
            end
         case 4912
            HEADER(img_indx).ANDOR.isBaselineClamped = entry.val;
         case 4913
            HEADER(img_indx).ANDOR.nPrescans = entry.val;
         case 4914
            HEADER(img_indx).ANDOR.model = entry.val;
         case 4915
            HEADER(img_indx).ANDOR.chipXSize = entry.val;
         case 4916
            HEADER(img_indx).ANDOR.chipYSize  = entry.val;
         case 4944
            HEADER(img_indx).ANDOR.baselineOffset = entry.val;
            
         case 33550       % GeoTIFF
            HEADER(img_indx).ModelPixelScaleTag = entry.val;
         case 33628       % Metamorph specific data
            HEADER(img_indx).MM_private1 = entry.val;
         case 33629       % this tag identify the image as a Metamorph stack!
            TIF.MM_stack = entry.val;
            TIF.MM_stackCnt = entry.cnt;
         case 33630       % Metamorph stack data: wavelength
            TIF.MM_wavelength = entry.val;
         case 33631       % Metamorph stack data: gain/background?
            TIF.MM_private2 = entry.val;
            
         case 33922       % GeoTIFF
            HEADER(img_indx).ModelTiePointTag = entry.val;
         case 34412       % Zeiss LSM data
            HEADER(img_indx).LSM_info = entry.val;
         case 34735       % GeoTIFF
            HEADER(img_indx).GeoKeyDirTag = entry.val;
         case 34737       % GeoTIFF
            HEADER(img_indx).GeoASCII = entry.val;
         case 42113       % GeoTIFF
            HEADER(img_indx).GDAL_NODATA = entry.val;
            
         case 50838       % Adobe
            HEADER(img_indx).meta_data_byte_counts = entry.val;
         case 50839       % Adobe
            HEADER(img_indx).meta_data = entry.val;
            
         otherwise
            if opt.ReadUnknownTags
               HEADER(img_indx).(['tag', num2str(entry_tag)])=entry.val;
               %eval(['IMG(img_indx+1).Tag',num2str(entry_tag),'=',entry.val,';']);
            else
               fprintf( 'Unknown TIFF entry with tag %i (cnt %i)\n', entry_tag, entry.cnt);
            end
      end
      
      % calculate bounding box  if you've got the stuff
      if isfield(HEADER, 'ModelPixelScaleTag') && isfield(HEADER, 'ModelTiePointTag') && isfield(HEADER, 'height')&& isfield(HEADER, 'width'),
         HEADER(img_indx).North=HEADER(img_indx).ModelTiePointTag(5)-HEADER(img_indx).ModelPixelScaleTag(2)*HEADER(img_indx).ModelTiePointTag(2);
         HEADER(img_indx).South=HEADER(img_indx).North-HEADER(img_indx).height*HEADER(img_indx).ModelPixelScaleTag(2);
         HEADER(img_indx).West=HEADER(img_indx).ModelTiePointTag(4)+HEADER(img_indx).ModelPixelScaleTag(1)*HEADER(img_indx).ModelTiePointTag(1);
         HEADER(img_indx).East=HEADER(img_indx).West+HEADER(img_indx).width*HEADER(img_indx).ModelPixelScaleTag(1);
      end
      
   end
   
   
   %total number of bytes per image:
   if TIF.PlanarConfiguration == 1
      TIF.BytesPerPlane = TIF.SamplesPerPixel * HEADER(img_indx).width * HEADER(img_indx).height * TIF.BytesPerSample;
   else
      TIF.BytesPerPlane = HEADER(img_indx).width * HEADER(img_indx).height * TIF.BytesPerSample;
   end
   
   % try to consolidate the TIFF strips if possible
   
   if opt.ConsolidateStrips
      consolidate_strips(TIF.BytesPerPlane);
   end
end

% determine the type of data stored in the pixels:
SampleFormat = 1;
if isfield(TIF, 'SampleFormat')
   SampleFormat = TIF.SampleFormat(1);
end

switch( SampleFormat )
   case 1
      TIF.classname = sprintf('uint%i', TIF.BitsPerSample(1));
   case 2
      TIF.classname = sprintf('int%i', TIF.BitsPerSample(1));
   case 3
      if (TIF.BitsPerSample(1) == 32 )
         TIF.classname = 'single';
      else
         TIF.classname = 'double';
      end
   otherwise
      error('unsuported TIFF sample format %i', SampleFormat);
end


%% Sub Functions

   function file_seek(fpos)
      status = fseek(TIF.file, fpos, -1);
      if status == -1
         error('Invalid file offset (invalid fseek)');
      end
   end


   function consolidate_strips(BytesPerPlane)
      
      for (nFrame = 1:numel(HEADER))
         
         %Try to consolidate the strips into a single one to speed-up reading:
         nBytes = HEADER(nFrame).StripByteCounts(1);
         
         if nBytes < BytesPerPlane
            
            idx = 1;
            % accumulate continguous strips that contribute to the plane
            while HEADER(nFrame).StripOffsets(1) + nBytes == HEADER(nFrame).StripOffsets(idx+1)
               idx = idx + 1;
               nBytes = nBytes + HEADER(nFrame).StripByteCounts(idx);
               if ( nBytes >= BytesPerPlane ); break; end
            end
            
            %Consolidate the Strips
            if ( nBytes <= BytesPerPlane(1) ) && ( idx > 1 )
               %fprintf('Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
               HEADER(nFrame).StripByteCounts = [nBytes; HEADER(nFrame).StripByteCounts(idx+1:HEADER(nFrame).StripNumber ) ];
               HEADER(nFrame).StripOffsets = HEADER(nFrame).StripOffsets( [1 , idx+1:HEADER(nFrame).StripNumber] );
               HEADER(nFrame).StripNumber  = 1 + HEADER(nFrame).StripNumber - idx;
            end
         end
      end
   end


%% ==================sub-functions that reads an IFD entry:===================


   function [nbBytes, matlabType] = convertType(tiffType)
      switch (tiffType)
         case 1 %byte
            nbBytes=1;
            matlabType='uint8';
         case 2 %ascii string
            nbBytes=1;
            matlabType='uchar';
         case 3 % word
            nbBytes=2;
            matlabType='uint16';
         case 4 %dword/uword
            nbBytes=4;
            matlabType='uint32';
         case 5 % rational
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
            error('tiff type %i not supported', tiffType)
      end
   end

%% ==================sub-functions that reads an IFD entry:===================

   function  entry = readIFDentry(entry_tag)
      
      entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
      entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      %disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);
      
      [ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);
      
      if entry.nbBytes * entry.cnt > 4
         %next field contains an offset:
         fpos = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
         %disp(strcat('offset = ', num2str(offset)));
         file_seek(fpos);
      end
      
      
      if entry_tag == 33629   % metamorph stack plane specifications
         entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.ByteOrder);
      elseif entry_tag == 34412  %TIF_CZ_LSMINFO
         entry.val = readLSMinfo;
      else
         if entry.tiffType == 5 % TIFF 'rational' type
            val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.ByteOrder);
            entry.val = val(1:2:length(val)) ./ val(2:2:length(val));
         else
            entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.ByteOrder);
         end
      end
      
      if ( entry.tiffType == 2 );
         entry.val = char(entry.val');
      end
      
   end


%% ==============partial-parse of LSM info:

   function R = readLSMinfo()
      
      % Read part of the LSM info table version 2
      % this provides only very partial information, since the offset indicate that
      % additional data is stored in the file
      
      R.MagicNumber          = sprintf('0x%09X',fread(TIF.file, 1, 'uint32', TIF.ByteOrder));
      S.StructureSize        = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.DimensionX           = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.DimensionY           = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.DimensionZ           = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.DimensionChannels    = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.DimensionTime        = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.IntensityDataType    = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.ThumbnailX           = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.ThumbnailY           = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.VoxelSizeX           = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      R.VoxelSizeY           = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      R.VoxelSizeZ           = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      R.OriginX              = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      R.OriginY              = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      R.OriginZ              = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      R.ScanType             = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
      R.SpectralScan         = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
      R.DataType             = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetVectorOverlay  = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetInputLut       = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetOutputLut      = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetChannelColors  = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      R.TimeInterval         = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
      S.OffsetChannelDataTypes = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetScanInformatio = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetKsData         = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetTimeStamps     = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetEventList      = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetRoi            = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetBleachRoi      = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      S.OffsetNextRecording  = fread(TIF.file, 1, 'uint32', TIF.ByteOrder);
      
      % There is more information stored in this table, which is skipped here
      
      %read real acquisition times:
      if ( S.OffsetTimeStamps > 0 )
         
         status =  fseek(TIF.file, S.OffsetTimeStamps, -1);
         if status == -1
            warning('tiffread:TimeStamps', 'Could not locate LSM TimeStamps');
            return;
         end
         
         StructureSize          = fread(TIF.file, 1, 'int32', TIF.ByteOrder);
         NumberTimeStamps       = fread(TIF.file, 1, 'int32', TIF.ByteOrder);
         for i=1:NumberTimeStamps
            R.TimeStamp(i)     = fread(TIF.file, 1, 'float64', TIF.ByteOrder);
         end
         
         %calculate elapsed time from first acquisition:
         R.TimeOffset = R.TimeStamp - R.TimeStamp(1);
         
      end
      
      % anything else assigned to S is discarded
      
   end

end


