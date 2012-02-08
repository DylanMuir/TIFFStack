% TIFFStack - Manipulate a TIFF file like a tensor
% 
% Usage: tsStack = TIFFStack(strFilename <, bInvert>)
% 
% A TIFFStack object behaves like a read-only memory mapped TIF file.  The
% entire image stack is treated as a matlab tensor.  Each frame of the file must
% have the same dimensions.  Reading the image data is optimised to the extent
% possible; the header information is only read once.
% 
% This class uses a modified version of tiffread [1, 2] to read data.  Code is
% included (but disabled) to use the matlab imread function, but this function
% returns invalid data for some TIFF formats.
% 
% Construction:
% 
% >> tsStack = TIFFStack('test.tiff');       % Construct a TIFF stack associated with a file
% 
% >> tsStack = TIFFStack('test.tiff', true); % Indicate that the image data should be inverted
% 
% tsStack = 
% 
%   TIFFStack handle
% 
%   Properties:
%          bInvert: 0
%      strFilename: [1x9 char]
%       sImageInfo: [5x1 struct]
%     strDataClass: 'uint16'
% 
% Usage:
% 
% >> tsStack(:, :, 3);     % Retrieve the 3rd frame of the stack, all planes
% 
% >> tsStack(:, :, 1, 3);  % Retrieve the 3rd plane of the 1st frame
% 
% >> size(tsStack)         % Find the size of the stack (rows, cols, frames, planes per pixel)
% 
% ans =
% 
%    128   128     5     1
% 
% >> tsStack(4);           % Linear indexing is supported
% 
% >> tsStack.bInvert = true;  % Turn on data inversion
% 
% References:
% [1] Francois Nedelec, Thomas Surrey and A.C. Maggs. Physical Review Letters
%        86: 3192-3195; 2001. DOI: 10.1103/PhysRevLett.86.3192
% 
% [2] http://www.embl.de/~nedelec/

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 28th June, 2011

classdef TIFFStack < handle
   properties
      bInvert;             % - A boolean flag that determines whether or not the image data will be inverted
   end
   
   properties (SetAccess = private)
      strFilename = [];    % - The name of the TIFF file on disk
      sImageInfo;          % - The TIFF header information
      strDataClass;        % - The matlab class in which data will be returned
   end
   
   properties (SetAccess = private, GetAccess = private)
      vnStackSize;
      TIF;                 % \_ Cached header infor for tiffread29 speedups
      HEADER;              % /
   end
   
   methods
      % TIFFStack - CONSTRUCTOR
      function oStack = TIFFStack(strFilename, bInvert)
         % - Check for inversion flag
         if (~exist('bInvert', 'var'))
            bInvert = false;
         end
         oStack.bInvert = bInvert;
         
         % - See if filename exists
         if (~exist(strFilename, 'file'))
            error('TIFFStack:InvalidFile', ...
                  '*** TIFFStack: File [%s] does not exist.', strFilename);
         end
         
         % - Assign absolute file path to stack
         oStack.strFilename = get_full_file_path(strFilename);
         
         % - Get image information
         try
            % - Read and save image information
            sInfo = imfinfo(strFilename);
            oStack.sImageInfo = sInfo;

            % - Read TIFF header for tiffread29
            [oStack.TIF, oStack.HEADER] = tiffread29_header(strFilename);
            
            % - Use imread to get the data class for this tiff
            % fPixel = imread(strFilename, 'TIFF', 1, 'PixelRegion', {[1 1], [1 1]});
            % oStack.strDataClass = class(fPixel);
            
            % - Use tiffread29 to get the data class for this tiff
            fPixel = tiffread29_readimage(oStack.TIF, oStack.HEADER, 1);
            fPixel = fPixel(1, 1, :);
            oStack.strDataClass = class(fPixel);

            % - Record stack size
            oStack.vnStackSize = [sInfo(1).Height sInfo(1).Width numel(sInfo) numel(fPixel)];
            
         catch mErr
            base_ME = MException('TIFFStack:InvalidFile', ...
                  '*** TIFFStack: Could not open file [%s].', strFilename);
            new_ME = addCause(base_ME, mErr);
            throw(new_ME);
         end
      end
      
      % delete - DESTRUCTOR
      function delete(oStack)
         % - Close the TIFF file, if opened by tiffread29_header
         if (isfield(oStack.TIF, 'file'))
            fclose(oStack.TIF.file);
         end
      end

%% --- Overloaded subsref

      function [tfData] = subsref(oStack, S)
         switch S(1).type
            case '()'
               nNumDims = numel(S.subs);
               nNumStackDims = numel(oStack.vnStackSize);
               
               % - Check dimensionality and trailing dimensions
               if (nNumDims == 1)
                  % - Translate from linear refs to indices
                  nNumDims = nNumStackDims;
                  
                  % - Translate colon indexing
                  if (isequal(S.subs{1}, ':'))
                     S.subs{1} = (1:prod(oStack.vnStackSize))';
                  end
                  
                  % - Get equivalent subscripted indexes
                  [S.subs{1:nNumDims}] = ind2sub(oStack.vnStackSize, S.subs{1});
                  
               elseif (nNumDims < nNumStackDims)
                  % - Assume trailing references are ':"
                  S.subs(nNumDims+1:nNumStackDims) = {':'};
                  
               elseif (nNumDims > nNumStackDims)
                  % - Check for non-colon references
                  vbNonColon = cellfun(@(c)(~ischar(c) | ~isequal(c, ':')), S.subs);
                  
                  % - Check only trailing dimensions
                  vbNonColon(1:nNumStackDims) = false;
                  
                  % - Check trailing dimensions for non-'1' indices
                  if (any(cellfun(@(c)(~isequal(c, 1)), S.subs(vbNonColon))))
                     % - This is an error
                     error('TIFFStack:badsubscript', ...
                        '*** TIFFStack: Index exceeds stack dimensions.');
                  end
                  
                  % - Only keep relevant dimensions
                  S.subs = S.subs(1:nNumStackDims);
               end
               
               % - Access stack
               tfData = TS_read_data_tiffread(oStack, S.subs);
               
               % - Reshape return data to concatenate trailing dimensions (just as
               % matlab does)
               if (nNumDims < nNumStackDims)
                  cnSize = num2cell(size(tfData));
                  tfData = reshape(tfData, cnSize{1:nNumDims-1}, []);
               end
               
            case '.'
               tfData = builtin('subsref', oStack, S);
               
            otherwise
               error('TIFFStack:InvalidReferencing', ...
                     '*** TIFFStack: Only ''()'' referencing is supported by TIFFStacks.');
         end
      end
      
%% --- Overloaded size
      function [varargout] = size(oStack, vnDimensions)
         % - Return the size of the stack
         vnSize = oStack.vnStackSize;
         
         % - Return specific dimension(s)
         if (exist('vnDimensions', 'var'))
            if (~isnumeric(vnDimensions))
               error('TIFFStack:dimensionMustBePositiveInteger', ...
                  '*** TIFFStack: Dimensions argument must be a positive integer within indexing range.');
            end
            
            % - Return the specified dimension(s)
            vnSize = vnSize(vnDimensions);
         end
         
         % - Handle differing number of size dimensions and number of output
         % arguments
         nNumArgout = max(1, nargout);
         
         if (nNumArgout == 1)
            % - Single return argument -- return entire size vector
            varargout{1} = vnSize;
            
         elseif (nNumArgout <= numel(vnSize))
            % - Several return arguments -- return single size vector elements,
            % with the remaining elements grouped in the last value
            varargout(1:nNumArgout-1) = num2cell(vnSize(1:nNumArgout-1));
            varargout{nNumArgout} = prod(vnSize(nNumArgout:end));
            
         else %(nNumArgout > numel(vnSize))
            % - Output all size elements
            varargout(1:numel(vnSize)) = num2cell(vnSize);
            
            % - Deal out trailing dimensions as '1'
            varargout(numel(vnSize)+1:nNumArgout) = {1};
         end
      end
      
%% --- Property accessors

      % set.bInvert - SETTER method for 'bInvert'
      function set.bInvert(oStack, bInvert)
         % - Check contents
         if (~islogical(bInvert) || ~isscalar(bInvert))
            error('TIFFStack:invalidArgument', ...
                  '*** TIFFStack/set.bInvert: ''bInvert'' must be a logical scalar.');
         else
            % - Assign bInvert value
            oStack.bInvert = bInvert;
         end
      end
      
   end
end

%% --- Helper functions ---

% TS_read_data_imread - FUNCTION Read the requested pixels from the TIFF file (using imread)
%
% Usage: [tfData] = TS_read_data_imread(oStack, cIndices)
%
% 'oStack' is a TIFFStack.  'cIndices' are the indices passed in from subsref.
% Colon indexing will be converted to full range indexing.  Reading is optimsed,
% by only reading pixels in a minimal-size window surrouding the requested
% pixels.  cIndices is a cell array with the format {rows, cols, frames,
% slices}.  Slices are RGB or CMYK or so on.

function [tfData] = TS_read_data_imread(oStack, cIndices) %#ok<DEFNU>
   % - Convert colon indexing
   vbIsColon = cellfun(@(c)(isequal(c, ':')), cIndices);
   
   for (nColonDim = find(vbIsColon))
      cIndices{nColonDim} = 1:oStack.vnStackSize(nColonDim);
   end
      
   % - Check ranges
   vnMinRange = cellfun(@(c)(min(c)), cIndices);
   vnMaxRange = cellfun(@(c)(max(c)), cIndices);
   
   if (any(vnMinRange < 1) || any(vnMaxRange > oStack.vnStackSize))
      error('TIFFStack:badsubscript', ...
            '*** TIFFStack: Index exceeds stack dimensions.');
   end
   
   % - Allocate large tensor
   vnBlockSize = vnMaxRange(1:2) - vnMinRange(1:2) + [1 1];
   vnBlockSize(3) = numel(cIndices{3});
   vnBlockSize(4) = oStack.vnStackSize(4);
   tfDataBlock = zeros(vnBlockSize, oStack.strDataClass);
   cBlock = {[vnMinRange(1) vnMaxRange(1)] [vnMinRange(2) vnMaxRange(2)]};
   
   % - Loop over frames to read data block
   try
      % - Disable warnings
      wOld = warning('off', 'MATLAB:rtifc:notPhotoTransformed');
      
      for (nFrame = 1:numel(cIndices{3})) %#ok<FORPF>
         tfDataBlock(:, :, nFrame, :) = imread(oStack.strFilename, 'TIFF', cIndices{3}(nFrame), 'PixelRegion', cBlock);
      end
      
      % - Re-enable warnings
      warning(wOld);
      
   catch mErr
      % - Record error state
      base_ME = MException('TIFFStack:ReadError', ...
                           '*** TIFFStack: Could not read data from image file.');
      new_ME = addCause(base_ME, mErr);
      throw(new_ME);
   end
   
   % - Do we need to resample the data block?
   bResample = any(~vbIsColon(1:3));
   if (bResample)
      cIndices{1} = cIndices{1} - vnMinRange(1) + 1;
      cIndices{2} = cIndices{2} - vnMinRange(2) + 1;
      tfData = tfDataBlock(cIndices{1}, cIndices{2}, :, cIndices{4});
   else
      tfData = tfDataBlock;
   end
   
   % - Invert data if requested
   if (oStack.bInvert)
      tfData = oStack.sImageInfo(1).MaxSampleValue - (tfData - oStack.sImageInfo(1).MinSampleValue);
   end
end


% TS_read_data_tiffread - FUNCTION Read the requested pixels from the TIFF file (using tiffread29)
%
% Usage: [tfData] = TS_read_data_imread(oStack, cIndices)
%
% 'oStack' is a TIFFStack.  'cIndices' are the indices passed in from subsref.
% Colon indexing will be converted to full range indexing.  cIndices is a cell
% array with the format {rows, cols, frames, slices}.  Slices are RGB or CMYK
% or so on.

function [tfData] = TS_read_data_tiffread(oStack, cIndices)
   % - Convert colon indexing
   vbIsColon = cellfun(@(c)(ischar(c) & isequal(c, ':')), cIndices);
   
   for (nColonDim = find(vbIsColon))
      cIndices{nColonDim} = 1:oStack.vnStackSize(nColonDim);
   end
      
   % - Check ranges
   vnMinRange = cellfun(@(c)(min(c)), cIndices);
   vnMaxRange = cellfun(@(c)(max(c)), cIndices);
   
   if (any(vnMinRange < 1) || any(vnMaxRange > oStack.vnStackSize))
      error('TIFFStack:badsubscript', ...
            '*** TIFFStack: Index exceeds stack dimensions.');
   end
   
   % - Allocate large tensor
   vnBlockSize = oStack.vnStackSize(1:2);
   vnBlockSize(3) = numel(cIndices{3});
   vnBlockSize(4) = oStack.vnStackSize(4);
   tfDataBlock = zeros(vnBlockSize, oStack.strDataClass);
   
   % - Read data block
   try
      tfDataBlock = tiffread29_readimage(oStack.TIF, oStack.HEADER, cIndices{3});
           
   catch mErr
      % - Record error state
      base_ME = MException('TIFFStack:ReadError', ...
                           '*** TIFFStack: Could not read data from image file.');
      new_ME = addCause(base_ME, mErr);
      throw(new_ME);
   end
   
   % - Do we need to resample the data block?
   bResample = any(~vbIsColon(1:3));
   if (bResample)
      tfData = tfDataBlock(cIndices{1}, cIndices{2}, :, cIndices{4});
   else
      tfData = tfDataBlock;
   end
   
   % - Invert data if requested
   if (oStack.bInvert)
      tfData = oStack.sImageInfo(1).MaxSampleValue - (tfData - oStack.sImageInfo(1).MinSampleValue);
   end
end

% get_full_file_path - FUNCTION Calculate the absolute path to a given (possibly relative) filename
%
% Usage: strFullPath = get_full_file_path(strFile)
%
% 'strFile' is a filename, which may include relative path elements.  The file
% does not have to exist.
%
% 'strFullPath' will be the absolute path to the file indicated in 'strFile'.

function strFullPath = get_full_file_path(strFile)

   [strDir, strName, strExt] = fileparts(strFile);

   if (isempty(strDir))
      strDir = '.';
   end

   strFullDirPath = cd(cd(strDir));
   strFullPath = fullfile(strFullDirPath, [strName strExt]);
end
