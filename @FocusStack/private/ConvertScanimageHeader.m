function sHeader = ConvertScanimageHeader(fileName)
% Dummy header conversion for FocusStack. Uses Scanimage internal functions
% to parse tif header information. Includes assignments2StructOrObj() and
% parseHeader() from Scanimage 3.81 distribution for ensured compatibility.
%
% Usage: [sHeader] = ConvertScanimageHeader(tif_filename)
%
% Ingie Hong <ingiehong@jhmi.edu>, 2015

%% Read TIFF file; extract # frames & image header
if ~exist(fileName,'file')
    error('''%s'' is not a recognized flag or filename. Aborting.',fileName);
end

warning('off','MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
hTif = Tiff(fileName);

headerString  = hTif.getTag('ImageDescription');

% numImages = 1;
% while ~hTif.lastDirectory()
%     numImages = numImages + 1;
%     hTif.nextDirectory();
% end
hTif.setDirectory(1);

if strncmp('state',headerString,5) 
    fileVersion = 3;
    header = parseHeader(headerString);
else
    fileVersion = 4;
    header = assignments2StructOrObj(headerString);            
end

%Extracts header info required by scim_openTif()
hdr = extractHeaderData(header, fileVersion);

%% Convert header information into FocusStack format

sHeader.vnFrameSizePixels = [hdr.numPixels hdr.numLines];
sHeader.tLineScanTime_ms = 1/ header.acq.frameRate / sHeader.vnFrameSizePixels(2) / 1e-3;
sHeader.vfXYZStep_nm = [0 0 header.acq.zStepSize*1000];
sHeader.fZoomFactor = [];  % - No stack resolution information is available
sHeader.nNumFrames = header.acq.numberOfZSlices * hdr.numFrames; 

sHeader.uNumChannels = hdr.numSlices * length(hdr.savedChans);

% - Visual stimulus information
if isfield(header, 'behavior') % if stimulus-related header information exists,
    sHeader.nStimulusID = header.behavior.nStimulusID;
    sHeader.tBlankTime = header.behavior.tBlankTime;
    %sHeader.vnSequenceIDs = str2num(header.behavior.StimSeqIDList) ; % gives error when empty
    eval(['sHeader.vnSequenceIDs = ' header.behavior.StimSeqIDList ' ; ' ]);
    sHeader.nSequenceLength = numel(sHeader.vnSequenceIDs);
else
    sHeader.vnSequenceIDs = [];
    sHeader.nSequenceLength = [];
    sHeader.nStimulusID = [];
    sHeader.tBlankTime = [];
end

% Copied from ConvertHelioscanHeader()
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
%sHeader.tBlankTime = [];
sHeader.nDataOffset = [];
sHeader.nDataLength = [];
sHeader.nScanLineLength = [];
sHeader.nZSinLength = [];

% - Order header fields
sHeader = orderfields(sHeader);

end


    function s = extractHeaderData(header, fileVersion)
       % Constants/Inits
        maxNumChans = 4;
        
        
        if fileVersion == 3
            localHdr = header;
            
            s.savedChans = [];
            for i=1:maxNumChans
                if isfield(localHdr.acq,['savingChannel' num2str(i)])
                    if localHdr.acq.(['savingChannel' num2str(i)]) && localHdr.acq.(['acquiringChannel' num2str(i)])
                        s.savedChans = [s.savedChans i];
                    end
                end
            end
            
            s.numPixels = localHdr.acq.pixelsPerLine;
            s.numLines = localHdr.acq.linesPerFrame;
            
            if isfield(localHdr.acq,'slowDimDiscardFlybackLine') && localHdr.acq.slowDimDiscardFlybackLine
                s.numLines = s.numLines - 1;
            end
            
            s.numSlices = localHdr.acq.numberOfZSlices;
            
            if ~localHdr.acq.averaging
                s.numFrames = localHdr.acq.numberOfFrames;
            else
                if isfield(localHdr.acq,'numAvgFramesSave')
                    s.numFrames = localHdr.acq.numberOfFrames / localHdr.acq.numAvgFramesSave;
                else
                    s.numFrames = 1;
                end
            end
            
            if  ~isfield(localHdr.internal,'lowPixelValue1')
                s.acqLUT = {};
            else
                s.acqLUT = cell(1,maxNumChans);
                for i=1:length(s.acqLUT)
                    s.acqLUT{i} = [localHdr.internal.(['lowPixelValue' num2str(i)]) localHdr.internal.(['highPixelValue' num2str(i)])];
                end
            end            
            
        elseif fileVersion == 4
            if isfield(header,'SI4App')
                localHdr = header.SI4App;
            else
                localHdr = header.SI4;
            end
            
            s.savedChans = localHdr.channelsSave;
            s.numPixels = localHdr.scanPixelsPerLine;
            s.numLines = localHdr.scanLinesPerFrame;
            s.acq.frameRate = localHdr.scanFrameRate;
            
            if isfield(localHdr,'acqNumAveragedFramesSaved')
                saveAverageFactor = localHdr.acqNumAveragedFramesSaved;
            elseif isfield(localHdr,'acqNumAveragedFrames')
                saveAverageFactor = localHdr.acqNumAveragedFrames;
            else
                assert(false);
            end

            s.numFrames = localHdr.acqNumFrames / saveAverageFactor;
            
            s.numSlices = localHdr.stackNumSlices;
            
            s.acqLUT = cell(1,size(localHdr.channelsLUT,1));
            for i=1:length(s.acqLUT)
                s.acqLUT{i} = localHdr.channelsLUT(i,:);
            end                
            
        else 
            assert(false);
        end

    end

    
function header = parseHeader(input)
% PARSEHEADER   - Read ScanImage Header String and return structure.
%   PARSEHEADER will output the value of the header fields as a structure.
%   Input is the header from ScanImage TIF File (char array).
%
% See also

out={};
tempcell=strread(input,'%q'); 
for lineCounter=1:length(tempcell)
    data=tempcell{lineCounter};
    if ~strncmp(data,'state.',6)
        out{end}=[out{end} ' ' data(1:end-1)];
        continue
    end
    equal=findstr('=',data);
    param=data(7:equal-1);
    val=data(equal+1:end);
    if isempty(val)
        val=[];
    elseif ~strcmp(val(1),'''')
        val=str2num(val);
    else
        if strcmp(val(end),'''')
            val=val(2:end-1);
        else
            val=val(2:end);
        end
    end
    out=[out {param} {val}];
end

while length(out)>2
    eval(['header.' out{1} '=out{2};']);
    out=out(3:end);
end
eval(['header.' out{1} '=out{2};']); % Added to catch last line.
end


function s = assignments2StructOrObj(str,s)
%ASSIGNMENTS2STRUCTOROBJ Create a struct, or configure a handle object,
%from a string generated by structOrObj2Assignments.
% s = assignments2Struct(str,s)
% 
% str: string generated by structOrObj2Assignments
% s [input]: (optional, scalar handle object). Object to configure.
% s [output]: If a handle object is supplied as s, that object is returned
% in s. If no object is supplied, a structure is created and returned in s.

if nargin < 2
    s = struct();
end

rows = textscan(str,'%s','Delimiter','\n');
rows = rows{1};

if isempty(rows)
    return;
end

for c = 1:numel(rows)
    row = rows{c};
    
    % replace top-level name with 'obj'
    [~, rmn] = strtok(row,'.');
    row = ['s' rmn];
    
    % deal with nonscalar nested structs/objs
    pat = '([\w]+)__([0123456789]+)\.';
    replc = '$1($2).';
    row = regexprep(row,pat,replc);
    
    % handle unencodeable value or nonscalar struct/obj.
    % Note: structOrObj2Assignments, assignments2StructOrObj, and toString
    % (all in most.util) are in cahoots with respect
    % to these hardcoded strings.
    unencodeval = '<unencodeable value>';
    if strfind(row,unencodeval)
        row = strrep(row,unencodeval,'[]');
    end
    nonscalarstructobjstr = '<nonscalar struct/object>';
    if strfind(row,nonscalarstructobjstr)
        row = strrep(row,nonscalarstructobjstr,'[]');
    end
    
    % handle ND array format produced by array2Str
    try 
        if ~isempty(strfind(row,'&'))
            equalsIdx = strfind(row,'=');
            [dimArr rmn] = strtok(row(equalsIdx+1:end),'&');
            arr = strtok(rmn,'&');
            arr = reshape(str2num(arr),str2num(dimArr)); %#ok<NASGU,ST2NM>
            eval([row(1:equalsIdx+1) 'arr;']);
        else
            eval([row ';']);
        end
    catch ME %Warn if assignments to no-longer-extant properties are found
        if strcmpi(ME.identifier,'MATLAB:noPublicFieldForClass')
            equalsIdx = strfind(row,'=');
            fprintf(1,'WARNING: Property ''%s'' was specified, but does not exist for class ''%s''\n', deblank(row(3:equalsIdx-1)),class(s));
        else
            ME.rethrow();
        end
    end
end

end

