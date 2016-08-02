function [TIF, INFO] = tiffread31_info(file_name)

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
INFO.SamplesPerPixel = 1;
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
   INFO.ByteOrder = 'little-endian';
elseif ( strcmp(bos','MM') )
   TIF.ByteOrder = 'ieee-be';
   INFO.ByteOrder = 'big-endian';
else
   error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies TIFF format

tiff_id = fread(TIF.file,1,'uint16', TIF.ByteOrder);

if (tiff_id ~= 42) && (tiff_id ~= 43)
   error('This is not a TIFF file (missing 42 or 43).');
end

% - By default, read 4-byte pointers
strIFDClassSize = 'uint32';
nIFDTagBytes = 12;
nIFDClassBytes = 2;
strTagSizeClass = 'uint32';
nInlineBytes = 4;

% - Handle a BigTIFF file
if (tiff_id == 43)
   offset_size = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
   test_val = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
   
   % - Test check value
   if (test_val ~= 0)
      error('This is not a valid BigTIFF file (invalid test value).');
   end
   
   % - Get IFD pointer data class
   switch (offset_size)
      case 8
         strIFDClassSize = 'uint64';
         nIFDClassBytes = 8;
         
      case 16
         strIFDClassSize = '2*uint64';
         nIFDClassBytes = 16;
         
      otherwise
         error('Unknown IFD pointer size for BigTIFF file.');
   end
   
   nIFDTagBytes = 20;
   strTagSizeClass = 'uint64';
   nInlineBytes = 8;
end

%% ---- read the image file directories (IFDs)

ifd_pos   = fread(TIF.file, 1, strIFDClassSize, TIF.ByteOrder);
img_indx  = 0;

while  ifd_pos ~= 0
   
   img_indx = img_indx + 1;
   
   TIF.file_name = file_name; %#ok<*AGROW>
   TIF.image_name = image_name;
   INFO(img_indx).index = img_indx;
   
   % move in the file to the next IFD
   file_seek(ifd_pos);
   
   %read in the number of IFD entries
   num_entries = fread(TIF.file,1,strIFDClassSize, TIF.ByteOrder);
   %fprintf('num_entries = %i\n', num_entries);
   
   % store current position:
   entry_pos = ifd_pos+nIFDClassBytes;
   
   % read the next IFD address:
   file_seek(ifd_pos+nIFDTagBytes*num_entries+nIFDClassBytes);
   ifd_pos = fread(TIF.file, 1, strIFDClassSize, TIF.ByteOrder);
   
   %fprintf('reading IFD %i at position %i\n', img_indx, ifd_pos);
   
   % read all the IFD entries
   for inx = 1:num_entries
      
      % move to next IFD entry in the file
      file_seek(entry_pos+nIFDTagBytes*(inx-1));
      
      % read entry tag
      entry_tag = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
      % read entry
      entry = readIFDentry(entry_tag, strTagSizeClass, nInlineBytes);
      
      switch entry_tag
         case 256         % image width = number of column
            INFO(img_indx).Width = entry.val;
         case 257         % image height = number of row
            INFO(img_indx).Height = entry.val;
            TIF.ImageLength = entry.val;
         case 258
            TIF.BitsPerSample = entry.val;
            TIF.BytesPerSample = TIF.BitsPerSample / 8;
            INFO(img_indx).BitsPerSample = entry.val';
         case 270         % general comments:
            INFO(img_indx).ImageDescription = entry.val;
         case 277
            TIF.SamplesPerPixel  = entry.val;
            INFO(img_indx).SamplesPerPixel  = entry.val;
            %fprintf('Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
         case 278
            INFO(img_indx).RowsPerStrip   = entry.val;
         case 284         %planar configuration describe the order of RGB
            TIF.PlanarConfiguration = entry.val;
         case 339
            TIF.SampleFormat   = entry.val;
      end
   end
   
   %total number of bytes per image:
   if TIF.PlanarConfiguration == 1
      TIF.BytesPerPlane = TIF.SamplesPerPixel * INFO(img_indx).Width * INFO(img_indx).Height * TIF.BytesPerSample;
   else
      TIF.BytesPerPlane = INFO(img_indx).Width * INFO(img_indx).Height * TIF.BytesPerSample;
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
      error('TIFFStack:Format', '*** TIFFStack: Error: Unsuported TIFF sample format [%i].', SampleFormat);
end


%% Sub Functions

   function file_seek(fpos)
      status = fseek(TIF.file, fpos, -1);
      if status == -1
         error('Invalid file offset (invalid fseek)');
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

         case 6 % signed byte
            nbBytes = 1;
            matlabType = 'int8';

         case 7   % undefined
            nbBytes=1;
            matlabType='uchar';
            
         case 8    % signed short
            nbBytes = 2;
            matlabType = 'int16';
            
         case 9   % signed long
            nbBytes = 4;
            matlabType = 'int32';
            
         case 10  % long rational
            nbBytes = 8;
            matlabType = 'int32';
         case 11
            nbBytes=4;
            matlabType='float32';
         case 12
            nbBytes=8;
            matlabType='float64';
            
         case 13 %TIFF_IFD
            nbBytes = 4;
            matlabType = 'uint32';
            
         case 16 % Long8
            nbBytes = 8;
            matlabType = 'uint64';
            
         case 17 % SLong8
            nbBytes = 8;
            matlabType = 'int64';
            
         case 18 % Unsigned IFD offet 8 bytes
            nbBytes = 8;
            matlabType = 'uint64';
            
         otherwise
            error('tiff type %i not supported', tiffType)
      end
   end

%% ==================sub-functions that reads an IFD entry:===================

   function  entry = readIFDentry(entry_tag, strTagSizeClass, nInlineBytes)
      
      entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.ByteOrder);
      entry.cnt      = fread(TIF.file, 1, strTagSizeClass, TIF.ByteOrder);
      %disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);
      
      [ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);
      
      if entry.nbBytes * entry.cnt > nInlineBytes
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
         
         StructureSize          = fread(TIF.file, 1, 'int32', TIF.ByteOrder); %#ok<NASGU>
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


