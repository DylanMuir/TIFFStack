function [sHeader] = ReadFocusHeader(hFileHandle)

% ReadFocusHeader - FUNCTION Read the header information from an FCS file
%
% Usage: [sHeader] = ReadFocusHeader(hFileHandle)
%
% 'hFileHandle' must be an open file handle.  Focus files must be opened in
% big-endian mode.
%


% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 9th December, 2010

% - Seek to file origin
fseek(hFileHandle, 0, 'bof');

% -- Read header
sHeader.fVersion = str2double(fread(hFileHandle, 1, 'uchar=>char')) + str2double(fread(hFileHandle, 1, 'uchar=>char')) / 10;
sHeader.strFilename = fread(hFileHandle, 14, 'uchar=>char')';
sHeader.UNUSED_strUsername = fread(hFileHandle, 16, 'uchar=>char')';
sHeader.strTime = fread(hFileHandle, 16, 'uchar=>char')';
sHeader.strDate = fread(hFileHandle, 16, 'uchar=>char')';

% - Trim fields
sHeader.strFilename = sHeader.strFilename(sHeader.strFilename ~= 0);
sHeader.UNUSED_strUsername = sHeader.UNUSED_strUsername(sHeader.UNUSED_strUsername ~= 0);
sHeader.strTime = sHeader.strTime(sHeader.strTime ~= 0);
sHeader.strDate = sHeader.strDate(sHeader.strDate ~= 0);

sHeader.tStartTime_ms = fread(hFileHandle, 1, 'uint32');
sHeader.tStopTime_ms = fread(hFileHandle, 1, 'uint32');

sHeader.vfXYZPositionStart_nm(3) = fread(hFileHandle, 1, 'int32');
sHeader.nHeaderLength = fread(hFileHandle, 1, 'uint32');
sHeader.vfXYZPositionStart_nm(1) = fread(hFileHandle, 1, 'int32');
sHeader.vfXYZPositionStart_nm(2) = fread(hFileHandle, 1, 'int32');
sHeader.vfXYZPositionStop_nm = fread(hFileHandle, 3, 'int32')';     % X Y Z

sHeader.UNUSED_fAttenuationDepth = fread(hFileHandle, 1, 'float32');

sHeader.bIntensityDevUsed = fread(hFileHandle, 1, 'uint8=>logical');
nul = fread(hFileHandle, 3, 'uchar');

sHeader.fStartIntensity = fread(hFileHandle, 1, 'float32');

sHeader.UNUSED_strReserved1 = fread(hFileHandle, 144, 'uchar=>char')';

sHeader.fBufferLength = fread(hFileHandle, 1, 'float32');
sHeader.UNUSED_fPixelBufferLength = fread(hFileHandle, 1, 'int32');
sHeader.UNUSED_vuFlags(1) = fread(hFileHandle, 1, 'uint32');
sHeader.UNUSED_uDAClockFreq = fread(hFileHandle, 1, 'uint32');

sHeader.fPixelClock_mHz = fread(hFileHandle, 1, 'float32');

sHeader.UNUSED_nXLeadTime_us = fread(hFileHandle, 1, 'int32');

sHeader.vnFrameSizePixels = fread(hFileHandle, 2, 'uint16')';

sHeader.UNUSED_vuFlags(2:3) = fread(hFileHandle, 2, 'uint16')';

sHeader.nPrePixels = fread(hFileHandle, 1, 'uchar');
sHeader.nPostPixels = fread(hFileHandle, 1, 'uchar');

sHeader.UNUSED_uYRetracePixels = fread(hFileHandle, 1, 'uint16');
sHeader.nNumFrames = round(fread(hFileHandle, 1, 'float32'));
sHeader.UNUSED_uNumChannels = fread(hFileHandle, 1, 'uint16');
sHeader.UNUSED_uSplitImages = fread(hFileHandle, 1, 'uint16');
sHeader.UNUSED_nYScanCutoff = fread(hFileHandle, 1, 'int32');

sHeader.UNUSED_uYLinesShutterOpenPreAcquisition = fread(hFileHandle, 1, 'uint16');
sHeader.UNUSED_tShutterOpenPreAquisition_us = fread(hFileHandle, 1, 'uint16');

sHeader.uNumAveragedFrames = fread(hFileHandle, 1, 'uint16');

sHeader.UNUSED_uToss = fread(hFileHandle, 1, 'uint16');

sHeader.vnFrameSizePixels_duplicate = fread(hFileHandle, 2, 'uint16')';
sHeader.tLineScanTime_ms = fread(hFileHandle, 1, 'float32');
sHeader.tRetraceDuration_ms = fread(hFileHandle, 1, 'float32');

sHeader.vfXYZStep_nm(3) = fread(hFileHandle, 1, 'int32');

sHeader.UNUSED_tWaitTimeZTSeriesScan = fread(hFileHandle, 1, 'uint32');

sHeader.UNUSED_vuFlags(4) = fread(hFileHandle, 1, 'uint32');
sHeader.UNUSED_strFlags = fread(hFileHandle, 4, 'uchar=>char')';
sHeader.UNUSED_vuFlags(5:6) = fread(hFileHandle, 2, 'int32')';

sHeader.vfXYZStep_nm(1:2) = fread(hFileHandle, 2, 'int32');

sHeader.UNUSED_strReserved2 = fread(hFileHandle, 156, 'uchar=>char')';

sHeader.vuXYScanRange = fread(hFileHandle, 2, 'uint16')';
sHeader.vuXYOffset = fread(hFileHandle, 2, 'uint16')';

sHeader.UNUSED_uParkingPosition = fread(hFileHandle, 1, 'uint16');

sHeader.uAngle_x100 = fread(hFileHandle, 1, 'uint16');

sHeader.UNUSED_vuFlags(7) = fread(hFileHandle, 1, 'int32');

sHeader.fZoomFactor = fread(hFileHandle, 1, 'float32');
sHeader.nNumFiles = round(fread(hFileHandle, 1, 'float32'));
sHeader.tTimeIntervalSec = fread(hFileHandle, 1, 'float32');
sHeader.fZRotation = fread(hFileHandle, 1, 'float32');

sHeader.UNUSED_strReserved3 = fread(hFileHandle, 24, 'uchar=>char')';

sHeader.UNUSED_uGain1 = fread(hFileHandle, 1, 'uint16');
sHeader.UNUSED_uOffset1 = fread(hFileHandle, 1, 'uint16');
sHeader.UNUSED_uGain2 = fread(hFileHandle, 1, 'uint16');
sHeader.UNUSED_uOffset2 = fread(hFileHandle, 1, 'uint16');

sHeader.vfXYPixelsFreeLine = fread(hFileHandle, 2, 'float32')';

sHeader.fPixelClockFreeLineHz = fread(hFileHandle, 1, 'float32');

sHeader.UNUSED_strReserved4 = fread(hFileHandle, 4, 'uchar=>char')';

sHeader.fAmplitudeFactor = fread(hFileHandle, 1, 'float32');
sHeader.fPixelClockZFreeLineHz = fread(hFileHandle, 1, 'float32');
sHeader.fPhaseShift = fread(hFileHandle, 1, 'float32');
sHeader.uModeSelector = fread(hFileHandle, 1, 'uint16');
sHeader.fCuspDelay = fread(hFileHandle, 1, 'float32');

sHeader.bFiberScope = fread(hFileHandle, 1, 'uchar=>logical');

sHeader.fPixelClockLSHz = fread(hFileHandle, 1, 'float32');
sHeader.fLSScansToRead = fread(hFileHandle, 1, 'float32');
sHeader.bLSXY = fread(hFileHandle, 1, 'uchar=>logical');

if (str2double(sHeader.strDate) < 101216)
   % - No StimServer information
   strReservedStimServer = fread(hFileHandle, 132, 'uchar=>char')';
   sHeader.nStimulusID = [];
   sHeader.tBlankTime = [];
   sHeader.nSequenceLength = [];
   sHeader.vnSequenceIDs = [];

elseif (str2double(sHeader.strDate) < 110106)
   % - Corrupted StimServer header
   sHeader.nStimulusID = fread(hFileHandle, 1, 'uint8');
   sHeader.tBlankTime = fread(hFileHandle, 1, 'single');
   sHeader.nSequenceLength = 80;
   sHeader.vnSequenceIDs = fread(hFileHandle, sHeader.nSequenceLength, 'uint8');
   nSkip = 132 - 1 - 4 - sHeader.nSequenceLength;
   fseek(hFileHandle, nSkip, 0);
   
   % - Fix up info
   sHeader.vnSequenceIDs = sHeader.vnSequenceIDs(sHeader.vnSequenceIDs ~= 0);
   sHeader.nSequenceLength = numel(sHeader.vnSequenceIDs);
   
else
   % - Correct StimServer header
   sHeader.nStimulusID = fread(hFileHandle, 1, 'uint8');
   sHeader.tBlankTime = fread(hFileHandle, 1, 'single');
   sHeader.nSequenceLength = fread(hFileHandle, 1, 'uint8');
   sHeader.vnSequenceIDs = fread(hFileHandle, sHeader.nSequenceLength, 'uint8');
   nSkip = 132 - 1 - 4 - 1 - sHeader.nSequenceLength;
   fseek(hFileHandle, nSkip, 0);
end

sHeader.nDataOffset = fread(hFileHandle, 1, 'int32');
sHeader.nDataLength = fread(hFileHandle, 1, 'int32');

sHeader.nScanLineLength = fread(hFileHandle, 1, 'int32');
sHeader.nZSinLength = fread(hFileHandle, 1, 'int32');

sHeader = orderfields(sHeader);


% -- Perform a basic sanity check on reading the header

if (any(sHeader.UNUSED_strReserved1)) || ...
      (any(sHeader.UNUSED_strReserved2)) || ...
      (any(sHeader.UNUSED_strReserved3)) || ...
      (any(sHeader.UNUSED_strReserved4)) || ...
      (sHeader.nDataOffset ~= 768)
   error('FocusStack:InvalidHeader', '*** FocusStack/ReadFocusHeader: Invalid header structure.');
end

% --- END of ReadFocusHeader.m ---

