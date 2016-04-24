% ConvertHelioscanHeader - PRIVATE FUNCTION Convert an XML header to Focus header structure format
%
% Usage: [sHeader] = ConvertHelioscanHeader(strHeaderXML)

function [sHeader] = ConvertHelioscanHeader(strHeaderXML)

% - Build a document model from the XML
SBISstream = java.io.StringBufferInputStream(strHeaderXML);
DBFfactory = javaMethod('newInstance', 'javax.xml.parsers.DocumentBuilderFactory');
DBbuilder = DBFfactory.newDocumentBuilder;
DOMheader = DBbuilder.parse(SBISstream);

% - Parse the XML to a structure
sHSHeader = parseNode(DOMheader);

% -- Extract header information

fDateNum = datenum(sHSHeader.OME.Image.AcquiredDate.Child_1.Data, 'yyyy-mm-ddTHH:MM:SS');
sHeader.strDate = datestr(fDateNum, 'yyyy-mm-dd');
sHeader.strTime = datestr(fDateNum, 'HH:MM:SS');

sPixels = sHSHeader.OME.Image.Pixels;

sHeader.vnFrameSizePixels = [str2double(sPixels.SizeX) str2double(sPixels.SizeY)];
sHeader.nNumFrames = str2double(sPixels.SizeT) * str2double(sPixels.SizeZ);

if (isfield(sHeader, 'tLineScanTime_ms'))
   sHeader.tLineScanTime_ms = (str2double(sPixels.TimeIncrement) * 1e3) / str2double(sPixels.SizeY);
else
   sHeader.tLineScanTime_ms = (1 / str2double(...
      sHSHeader.OME.Image.StructuredAnnotations.XmlStringAnnotation.HelioScan.ImagingMode.settings.FrameSettings.FrameScanRate.Child_1.Data) ...
      * 1e3) / str2double(sPixels.SizeY);
end      

sHeader.vfXYZStep_nm = [0 0 0];

sHeader.fZoomFactor = 117 / (str2double(sPixels.PhysicalSizeX));

if (isfield(sHSHeader.OME.Image.StructuredAnnotations.XmlStringAnnotation.HelioScan.Stimulator, 'settings'))
   strTCPSignalPath = sHSHeader.OME.Image.StructuredAnnotations.XmlStringAnnotation.HelioScan.Stimulator.settings.TCP_signal.Child_1.Data;
   sHeader.vnSequenceIDs = eval(strTCPSignalPath(strfind(strTCPSignalPath,'['):strfind(strTCPSignalPath,']')));
   sHeader.nSequenceLength = numel(sHeader.vnSequenceIDs);
   
   sHeader.nStimulusID = sscanf(strTCPSignalPath, 'SS PRESENT %d');
else
   sHeader.vnSequenceIDs = [];
   sHeader.nSequenceLength = [];
   sHeader.nStimulusID = [];
end

sHeader.fVersion = [];
sHeader.strFilename = [];
sHeader.tStartTime_ms = [];
sHeader.tStopTime_ms = [];
sHeader.vfXYZPositionStart_nm = [];    % X Y Z
sHeader.vfXYZPositionStop_nm = []';     % X Y Z
sHeader.nHeaderLength = [];
sHeader.UNUSED_fAttenuationDepth = [];
sHeader.bIntensityDevUsed = [];
sHeader.fStartIntensity = [];
sHeader.UNUSED_strReserved1 = [];
sHeader.fBufferLength = [];
sHeader.UNUSED_fPixelBufferLength = [];
sHeader.UNUSED_uDAClockFreq = [];
sHeader.fPixelClock_mHz = [];
sHeader.UNUSED_nXLeadTime_us = [];
sHeader.nPrePixels = [];
sHeader.nPostPixels = [];
sHeader.UNUSED_uYRetracePixels = [];
sHeader.UNUSED_uNumChannels = [];
sHeader.UNUSED_uSplitImages = [];
sHeader.UNUSED_nYScanCutoff = [];
sHeader.UNUSED_uYLinesShutterOpenPreAcquisition = [];
sHeader.UNUSED_tShutterOpenPreAquisition_us = [];
sHeader.uNumAveragedFrames = [];
sHeader.UNUSED_uToss = [];
sHeader.vnFrameSizePixels_duplicate = [];
sHeader.tRetraceDuration_ms = [];
sHeader.UNUSED_tWaitTimeZTSeriesScan = [];
sHeader.UNUSED_vuFlags = [];
sHeader.UNUSED_strFlags = [];
sHeader.UNUSED_strReserved2 = [];
sHeader.vuXYScanRange = [];
sHeader.vuXYOffset = [];
sHeader.UNUSED_uParkingPosition = [];
sHeader.uAngle_x100 = [];
sHeader.nNumFiles = [];
sHeader.tTimeIntervalSec = [];
sHeader.fZRotation = [];
sHeader.UNUSED_strReserved3 = [];
sHeader.UNUSED_uGain1 = [];
sHeader.UNUSED_uOffset1 = [];
sHeader.UNUSED_uGain2 = [];
sHeader.UNUSED_uOffset2 = [];
sHeader.vfXYPixelsFreeLine = [];
sHeader.fPixelClockFreeLineHz = [];
sHeader.UNUSED_strReserved4 = [];
sHeader.fAmplitudeFactor = [];
sHeader.fPixelClockZFreeLineHz = [];
sHeader.fPhaseShift = [];
sHeader.uModeSelector = [];
sHeader.fCuspDelay = [];
sHeader.bFiberScope = [];
sHeader.fPixelClockLSHz = [];
sHeader.fLSScansToRead = [];
sHeader.bLSXY = [];
sHeader.tBlankTime = [];
sHeader.nDataOffset = [];
sHeader.nDataLength = [];
sHeader.nScanLineLength = [];
sHeader.nZSinLength = [];


% - Order header fields
sHeader = orderfields(sHeader);

end


function sNode = parseNode(theNode)
% Recurse over node children.
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   
   for count = 1:numChildNodes
      % - Recursively parse child nodes
      sChild = parseNode(childNodes.item(count-1));
      
      % - Assign child node to this node
      try
         sNode.(char(childNodes.item(count-1).getNodeName)) = sChild;
         
      catch mErr
         % - Was it an invalid field name?
         if (isequal(mErr.identifier, 'MATLAB:AddField:InvalidFieldName'))
            % - Fix this by assigning the name within the child
            sChild.Name = char(childNodes.item(count-1).getNodeName);
            sNode.(sprintf('Child_%d', count)) = sChild;
            
         else
            % - Failed, so retrhow the error
            base_ME = mException('FocusStack:FailedConvertingHeliscanHeader', ...
               '*** FocusStack/ConvertHelioscanHeader: Could not parse XML header structure.');
            new_ME = addCause(base_ME, mErr);
            throw(new_ME);
         end
      end
   end
end

% - Extract attributes for this node
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   
   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      
      try
         sNode.(char(attrib.getName)) = char(attrib.getValue);
         
      catch mErr
         % - Was it an invalid field name?
         if (isequal(mErr.identifier, 'MATLAB:AddField:InvalidFieldName'))
            % - Fix this by assigning the name and value within the attribute
            sAttribute.Name = char(attrib.getName);
            sAttribute.Value = char(attrib.getValue);
            sNode.(sprintf('Attribute_%d', count)) = sAttribute;
            
         else
            % - Failed, so retrhow the error
            base_ME = mException('FocusStack:FailedConvertingHeliscanHeader', ...
               '*** FocusStack/ConvertHelioscanHeader: Could not parse XML header structure.');
            new_ME = addCause(base_ME, mErr);
            throw(new_ME);
         end
      end
   end
end

% - Extract data for this node
if any(strcmp(methods(theNode), 'getData'))
   sNode.Data = char(theNode.getData);
else
   sNode.Data = '';
end

end

