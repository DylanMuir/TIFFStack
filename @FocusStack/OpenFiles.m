function OpenFiles(oStack)

   % OpenFiles - METHOD Open the files for this stack

   nNumFiles = numel(oStack.cstrFilenames);

   try
      for (nFile = 1:nNumFiles)
         % - Extract a full file path
         strFullPath = get_full_file_path(oStack.cstrFilenames{nFile});
         [~, strFilenameOnly, strExt] = fileparts(strFullPath);
         strFilenameOnly = [strFilenameOnly strExt]; %#ok<AGROW>

         % - Check if file exists
         if (~exist(strFullPath, 'file'))
            error('FocusStack:FileNotFound', ...
               '*** FocusStack/OpenFiles: Could not find the specifed raw file [%s].', ...
               strFullPath);
         end

         % - Record full file path
         oStack.cstrFilenames{nFile} = strFullPath;

         % - Check raw data type
         switch lower(strExt)
            case {'.fcs'}
               OpenFocusStack(oStack, strFullPath, strFilenameOnly, nFile);

            case {'.tif', '.tiff'}
                % Check tif file header for file format
                % (Scanimage/Prarie etc)
                info=imfinfo(oStack.cstrFilenames{nFile});
                headerString = info(1).ImageDescription;
                if strncmp('state', headerString, 5) % If Scanimage 3.X format... this should be later updated to accomodate later versions and Prarie etc.
                    OpenScanimageTifStack(oStack, strFullPath, strFilenameOnly, nFile);
                elseif strfind(headerString, 'scanimage.SI4')
                    OpenScanimageTifStack(oStack, strFullPath, strFilenameOnly, nFile);
                else
                    OpenTifStack(oStack, strFullPath, strFilenameOnly, nFile);
                end
                
            case {'.bin'}
               OpenBinStack(oStack, strFullPath, strFilenameOnly, nFile);
         end

      end

   catch sErr
      % - Display an error message
      disp('*** FocusStack: Could not open the specified raw data files.');

      rethrow(sErr);
   end

end

function strFullPath = get_full_file_path(strFile)

   [strDir, strName, strExt] = fileparts(strFile);

   if (isempty(strDir))
      strDir = '.';
   end

   strFullDirPath = cd(cd(strDir));
   strFullPath = fullfile(strFullDirPath, [strName strExt]);
end

function OpenFocusStack(oStack, strFullFilePath, strFilenameOnly, nFile) %#ok<INUSL>
   % - Try to open the file (big-endian mode)
   [hFileHandle, strMessage] = fopen(oStack.cstrFilenames{nFile}, 'r', 'b');

   % - Could we open the file?
   if (hFileHandle == -1)
      error('FocusStack:FileNotOpened', ...
         '*** FocusStack/OpenFiles/OpenFocusStack: Could not open the specified focus raw file [%s]\n   Error message: %s.', ...
         ['.../' strFilenameOnly], strMessage);
   end

   % - Read focus file header and record it
   sHeader = ReadFocusHeader(hFileHandle);

   if (nFile == 1)
      oStack.vsHeaders = sHeader;
   else
      oStack.vsHeaders(nFile) = sHeader;
   end

   % - Check that this stack is compatible with previous stack sizes (frame
   % sizes)
   if (isempty(oStack.vnFrameSize))
      oStack.vnFrameSize = sHeader.vnFrameSizePixels;

   elseif (~isequal(sHeader.vnFrameSizePixels, oStack.vnFrameSize))
      error('FocusStack:DifferentFrameSizes', ...
         '*** FocusStack/OpenFiles/OpenFocusStack: Raw file [%s] has a different frame size than the stack.', ...
         ['.../' strFilenameOnly]);
   end

   % - Check frame duration
   tThisFrameDuration = sHeader.tLineScanTime_ms*1e-3 * sHeader.vnFrameSizePixels(2);

   if (isempty(oStack.tFrameDuration))
      oStack.tFrameDuration = tThisFrameDuration;

   elseif (~isequal(round(oStack.tFrameDuration * 1000), round(tThisFrameDuration * 1000)))
      warning('FocusStack:DifferentFrameDuration', ...
         '--- FocusStack/OpenFiles/OpenFocusStack: Raw file [%s] has a different frame duration than the stack (%dms vs %dms).', ...
         ['.../' strFilenameOnly], round(oStack.tFrameDuration * 1000), round(tThisFrameDuration * 1000));
   end

   % - Check Z step
   if (isempty(oStack.fZStep))
      oStack.fZStep = sHeader.vfXYZStep_nm(3) / 1e9;

   elseif (~isequal(oStack.fZStep, sHeader.vfXYZStep_nm(3) / 1e9))
      warning('FocusStack:DifferentZStep', ...
         '--- FocusStack/OpenFiles/OpenFocusStack: Raw file [%s] has a different frame Z step than the stack.', ...
         ['.../' strFilenameOnly]);
   end

   % - Check zoom
   if (isempty(oStack.fPixelsPerUM))
      oStack.fPixelsPerUM = sHeader.vnFrameSizePixels(1) / (117 / sHeader.fZoomFactor);

   elseif (~isequal(oStack.fPixelsPerUM, sHeader.vnFrameSizePixels(1) / (117 / sHeader.fZoomFactor)))
      warning('FocusStack:DifferentZoom', ...
         '--- FocusStack/OpenFiles/OpenFocusStack: Raw file [%s] has a different zoom level than the stack (%dum vs %dum).', ...
         ['.../' strFilenameOnly], round(oStack.vnFrameSize(1) / oStack.fPixelsPerUM), ...
         round(oStack.vnFrameSize(1) / (sHeader.vnFrameSizePixels(1) / (117 / sHeader.fZoomFactor))));
   end

   % - Set number of frames
   oStack.vnNumFrames(nFile) = sHeader.nNumFrames;

   % - Close file
   fclose(hFileHandle);

   % - Memory-map the file
%    oStack.vhMemMapFileHandles{nFile} = ...
%       memmapfile(oStack.cstrFilenames{nFile}, 'Format', {'uint8' [2 prod(oStack.vnFrameSize) oStack.vnNumFrames(nFile)] 'tfStack'}, 'Offset', 768);
   oStack.vhMemMapFileHandles{nFile} = MappedTensor(oStack.cstrFilenames{nFile}, 2, prod(oStack.vnFrameSize), oStack.vnNumFrames(nFile), 'Class', 'uint8', 'HeaderBytes', 768);

   % - Two channels for focus files
   oStack.nNumChannels = 2;
   
   % - Focus files are always uint8
   oStack.strDataClass = 'uint8';
end

function OpenTifStack(oStack, strFullPath, strFilenameOnly, nFile)
   % - Construct TIFFStack
   oStack.vhMemMapFileHandles{nFile} = TIFFStack(strFullPath);
   
   % - Transpose the stack to match old usage
   oStack.vhMemMapFileHandles{nFile} = transpose(oStack.vhMemMapFileHandles{nFile});
   
   % - Try to convert a Helioscan header, if it exists
   try
      sImageInfo = getImageInfo(oStack.vhMemMapFileHandles{nFile});
      sHeader = ConvertHelioscanHeader(sImageInfo(1).ImageDescription);
      
   catch
      % - Try to extract basic stack information
      vnThisStackSize = size(oStack.vhMemMapFileHandles{nFile});
      sHeader.vnFrameSizePixels = vnThisStackSize(1:2);
      sHeader.tLineScanTime_ms = oStack.tFrameDuration / sHeader.vnFrameSizePixels(2) / 1e-3;
      sHeader.vfXYZStep_nm(3) = 0;
      sHeader.fZoomFactor = 1./((oStack.fPixelsPerUM ./ sHeader.vnFrameSizePixels(1)) ./ 117);
      sHeader.nNumFrames = vnThisStackSize(3);
      sHeader.nStimulusID = nan;
      sHeader.tBlankTime = nan;

      % - Sequence is unknown
      sHeader.vnSequenceIDs = nan;
   end
   
   if (nFile == 1)
      oStack.vsHeaders = sHeader;
   else
      oStack.vsHeaders(nFile) = sHeader;
   end

   % - Check that this stack is compatible with previous stack sizes (frame
   % sizes)
   if (isempty(oStack.vnFrameSize))
      oStack.vnFrameSize = sHeader.vnFrameSizePixels;

   elseif (~isequal(sHeader.vnFrameSizePixels, oStack.vnFrameSize))
      error('FocusStack:DifferentFrameSizes', ...
         '*** FocusStack/OpenFiles/OpenTifStack: Raw file [%s] has a different frame size than the stack.', ...
         ['.../' strFilenameOnly]);
   end

   % - Check frame duration
   tThisFrameDuration = sHeader.tLineScanTime_ms*1e-3 * sHeader.vnFrameSizePixels(2);

   if (isempty(oStack.tFrameDuration))
      oStack.tFrameDuration = tThisFrameDuration;

   elseif (~isequal(round(oStack.tFrameDuration * 1000), round(tThisFrameDuration * 1000)))
      warning('FocusStack:DifferentFrameDuration', ...
         '--- FocusStack/OpenFiles/OpenTifStack: Raw file [%s] has a different frame duration than the stack (%dms vs %dms).', ...
         ['.../' strFilenameOnly], round(oStack.tFrameDuration * 1000), round(tThisFrameDuration * 1000));
   end

   % - Check Z step
   if (isempty(oStack.fZStep))
      oStack.fZStep = sHeader.vfXYZStep_nm(3) / 1e9;

   elseif (~isequal(oStack.fZStep, sHeader.vfXYZStep_nm(3) / 1e9))
      warning('FocusStack:DifferentZStep', ...
         '--- FocusStack/OpenFiles/OpenTifStack: Raw file [%s] has a different frame Z step than the stack.', ...
         ['.../' strFilenameOnly]);
   end

   % - Check zoom
   if (isempty(oStack.fPixelsPerUM))
      oStack.fPixelsPerUM = sHeader.vnFrameSizePixels(1) ./ (117 ./ sHeader.fZoomFactor);

   elseif (~isequal(oStack.fPixelsPerUM, sHeader.vnFrameSizePixels(1) ./ (117 ./ sHeader.fZoomFactor)))
      warning('FocusStack:DifferentZoom', ...
         '--- FocusStack/OpenFiles/OpenTifStack: Raw file [%s] has a different zoom level than the stack (%dum vs %dum).', ...
         ['.../' strFilenameOnly], round(oStack.vnFrameSize(1) ./ oStack.fPixelsPerUM), ...
         round(oStack.vnFrameSize(1) ./ (sHeader.vnFrameSizePixels(1) ./ (117 ./ sHeader.fZoomFactor))));
   end
   
   % - Check number of channels
   if (isempty(oStack.nNumChannels))
      oStack.nNumChannels = size(oStack.vhMemMapFileHandles{nFile}, 4);
      
   elseif (~isequal(oStack.nNumChannels, size(oStack.vhMemMapFileHandles{nFile}, 4)))
      error('FocusStack:DifferentFrameSizes', ...
            '*** FocusStack/OpenFiles/OpenTifStack: Raw file [%s] has a different number of channels than the stack.', ...
            ['.../' strFilenameOnly]);
   end

   % - Check data class
   if (isempty(oStack.strDataClass))
      oStack.strDataClass = getDataClass(oStack.vhMemMapFileHandles{nFile});
      
   elseif (~isequal(oStack.strDataClass, getDataClass(oStack.vhMemMapFileHandles{nFile})))
      error('FocusStack:DifferentDataClass', ...
            '*** FocusStack/OpenFiles/OpenTifStack: Raw file [%s] has a different data class than the stack.', ...
            ['.../' strFilenameOnly]);
   end
   
   % - Set number of frames for this file
   oStack.vnNumFrames(nFile) = sHeader.nNumFrames;
end

function OpenBinStack(oStack, strFullPath, strFilenameOnly, nFile)
   % - find header name
   [strDir, strName, ~] = fileparts(strFullPath);
   strHeaderName = strName(strfind(strName,'2') : strfind(strName,'h'));
   
   % - Try to convert a Helioscan header, if it exists
%    try
      sHeader = ConvertXMLHeader(fileread(fullfile(strDir, [strHeaderName '.xml'])));
      
%    catch mErr
%       error('FocusStack:ReadingXMLHeader', ...
%             '*** FocusStack/OpenFiles/OpenBinStack: XML header can not be read.', ...
%             ['.../' strFilenameOnly]);
%    end
   
   if (nFile == 1)
      oStack.vsHeaders = sHeader;
   else
      oStack.vsHeaders(nFile) = sHeader;
   end

   % - Construct BinStack
%    oStack.vhMemMapFileHandles{nFile} = MappedTensor(oStack.cstrFilenames{nFile}, prod(oStack.vsHeaders.vnFrameSizePixels), oStack.vsHeaders.nNumFrames, 'Class', 'uint16');
   oStack.vhMemMapFileHandles{nFile} = ...
      MappedTensor(  oStack.cstrFilenames{nFile}, [1 prod(oStack.vsHeaders(nFile).vnFrameSizePixels) oStack.vsHeaders(nFile).nNumFrames], ...
                     'Class', 'uint16', 'MachineFormat', 'ieee-be.l64');
%    oStack.vhMemMapFileHandles{nFile} = MappedTensor(oStack.cstrFilenames{nFile}, 1, prod(oStack.vsHeaders.vnFrameSizePixels), oStack.vsHeaders.nNumFrames, 'Class', 'uint16');

   % - Check that this stack is compatible with previous stack sizes (frame
   % sizes)
   if (isempty(oStack.vnFrameSize))
      oStack.vnFrameSize = sHeader.vnFrameSizePixels;

   elseif (~isequal(sHeader.vnFrameSizePixels, oStack.vnFrameSize))
      error('FocusStack:DifferentFrameSizes', ...
         '*** FocusStack/OpenFiles/OpenBinStack: Raw file [%s] has a different frame size than the stack.', ...
         ['.../' strFilenameOnly]);
   end

   % - Check frame duration
   tThisFrameDuration = sHeader.tFrameDuration;

   if (isempty(oStack.tFrameDuration))
      oStack.tFrameDuration = tThisFrameDuration;

   elseif (~isequal(round(oStack.tFrameDuration * 1000), round(tThisFrameDuration * 1000)))
      warning('FocusStack:DifferentFrameDuration', ...
         '--- FocusStack/OpenFiles/OpenBinStack: Raw file [%s] has a different frame duration than the stack (%dms vs %dms).', ...
         ['.../' strFilenameOnly], round(oStack.tFrameDuration * 1000), round(tThisFrameDuration * 1000));
   end

   % - Check Z step
   if (isempty(oStack.fZStep))
      oStack.fZStep = sHeader.vfXYZStep_nm(3) / 1e9;

   elseif (~isequal(oStack.fZStep, sHeader.vfXYZStep_nm(3) / 1e9))
      warning('FocusStack:DifferentZStep', ...
         '--- FocusStack/OpenFiles/OpenBinStack: Raw file [%s] has a different frame Z step than the stack.', ...
         ['.../' strFilenameOnly]);
   end

   % - Check zoom
   if (isempty(oStack.fPixelsPerUM))
      oStack.fPixelsPerUM = sqrt(prod(sHeader.vnFrameSizePixels ./ (sHeader.vfFieldOfView * 1e6))) * sHeader.fZoomFactor;

   elseif (~isequal(oStack.fPixelsPerUM, sqrt(prod(sHeader.vnFrameSizePixels ./ (sHeader.vfFieldOfView * 1e6))) * sHeader.fZoomFactor))
      warning('FocusStack:DifferentZoom', ...
         '--- FocusStack/OpenFiles/OpenBinStack: Raw file [%s] has a different zoom level than the stack (%dum vs %dum).', ...
         ['.../' strFilenameOnly], round(oStack.vnFrameSize(1) ./ oStack.fPixelsPerUM), ...
         round(oStack.vnFrameSize(1) / sqrt(prod(sHeader.vnFrameSizePixels ./ (sHeader.vfFieldOfView * 1e6))) * sHeader.fZoomFactor));
   end
   
   % - Check number of channels
   if (isempty(oStack.nNumChannels))
      oStack.nNumChannels = size(oStack.vhMemMapFileHandles{nFile}, 1);
      
   elseif (~isequal(oStack.nNumChannels, size(oStack.vhMemMapFileHandles{nFile}, 1)))
      error('FocusStack:DifferentFrameSizes', ...
            '*** FocusStack/OpenFiles/OpenBinStack: Raw file [%s] has a different number of channels than the stack.', ...
            ['.../' strFilenameOnly]);
   end

   % - Bin files are always uint16
   oStack.strDataClass = 'uint16';
   
   % - Set number of frames for this file
   oStack.vnNumFrames(nFile) = sHeader.nNumFrames;
end

% --- END of OpenFiles METHOD ---
