% MappedTensor - CLASS Create and map a new file to a large variable
%
% The MappedTensor class creates a large variable, and maps it directly to a
% file on disk.  The variable can be passed around BY REFERENCE, indexed and
% written to without allocating space for the entire variable in matlab.  Note
% that this is a handle class, meaning that if you copy the variable you only
% copy the handle and not the data.  Modifying one tensor will modify all the
% copies!
%
% Creation: mtVariable = MappedTensor(vnTensorSize)
%           mtVariable = MappedTensor(nDim1, nDim2, nDim3, ...)
%           mtVariable = MappedTensor(strExistingFilename, ...)
%           mtVariable = MappedTensor(..., 'Class', strClassName)
%           mtVariable = MappedTensor(..., 'HeaderBytes', nHeaderBytesToSkip)
%
% 'vnTensorSize', or [nDim1 nDim2 nDim3 ...] defines the desired size of the
% variable.  By default, a new binary temporary file will be generated, and
% deleted when the 'mtVariable' is destroyed.  'strExistingFilename' can be used
% to map an existing file on disk, but the full size (and class) of the file
% must be known and specified in advance.  This file will not be removed when
% all handle references are destroyed.
%
% By default the tensor will have class 'double'.  This can be specified as an
% argument to MappedTensor.  Supported classes: char, int8, uint8, logical,
% int16, uint16, int32, uint32, single, int64, uint64, double.
%
% The optional parameter 'nHeaderBytesToSkip' allows you to skip over the
% beginning of an (existing) binary file, by throwing away the specified
% number of header bytes.
%
% Usage: size(mtVariable)
%        mtVariable(:) = rand(100, 100, 100);
%        mfData = mtVariable(:, :, 34, 2);
%        mtVariable(12) = 1+6i;
%
% Note: mtVariable = rand(100, 100, 100); would over-write the mapped tensor
% with a standard matlab tensor!  To assign to the entire tensor you must use
% colon referencing: mtVariable(:) = ...
%
% It's not clear why you would do this anyway, because the right hand side of
% the assignment would already allocate enough space for the full tensor...
% which is presumably what you're trying to avoid.
%
% Permute is supported.  Complex numbers are supported (a definite benefit over
% memmapfile).  Transpose (.') and ctranspose (') are both supported.
% Transposition just swaps the first two dimensions, leaving the trailing
% dimensions unpermuted.
%
% Unary plus (+A) and minus (-A) are supported.  Binary plus (A+B), minus (A-B),
% times (A*B, A.*B) as long as one of A or B is a scalar.  Divide (A/B,
% A./B, B\A, B.\A) is supported, as long as B is a scalar.
%
% Save and load is minimally supported -- data is NOT saved, but on load a new
% mapped tensor will be generated and filled with zeros.  Both save and load
% generate warnings.
%
% Dot referencing ('.') is not supported.
%
% sum(mtVar <, nDimension>) is supported, to avoid importing the entire tensor
% into memory.
%
% Convenience functions:
%    SliceFunction: Execute a function on the entire tensor, by slicing it along
%       a specified dimension, and store the results back in the tensor.  
%
%       Usage: [<mtNewVar>] = SliceFunction(mtVar, fhFunctionHandle, nSliceDim <, vnSliceSize,> ...)
%
%           'mtVar' is a MappedTensor.  This tensor will be sliced up along
%           dimensions 'nSliceDim', with each slice passed individually to
%           'fhFunctionHandle', along with any trailing argments (...).  If no
%           return argument is supplied, the results will be stored back in
%           'mtVar'.  If a return argument is supplied, a new MappedTensor will
%           be created to contain the results.  The optional argument
%           'vnSliceSize' can be used to call a function that returns a
%           different sized output than the size of a single slice of 'mtVar'.
%           In that case, a new tensor 'mtNewVar' will be generated, and it will
%           have the size 'vnSliceSize', with the dimension 'nSliceDim' having
%           the same length as in the original tensor 'mtVar'.
%
%           "Slice assign" operations can be performed by passing in a function
%           than takes no input arguments for 'fhFunctionHandle'.
%
%       For example:
%
%       mtVar(:) = abs(fft2(mtVar(:, :, :)));
%
%          is equivalent to
%
%       SliceFunction(mtVar, @(x)(abs(fft2(x)), 3);
%
%       Each slice of the third dimension of mtVar, taken in turn, is passed to
%       fft2 and the result stored back into the same slice of mtVar.
%
%       mtVar2 = SliceFunction(mtVar, @fft2, 3);
%
%       This will return the result in a new MappedTensor, with temporary
%       storage.
%
%       mtVar2 = SliceFunction(mtVar, @sum, 3, [1 10 1]);
%
%       This will create a new MappedTensor with size [1 10 N], where 'N' is
%       the length along dimension 3 of 'mtVar'.
%
%       SliceFunction(mtVar, @()(randn(10, 10)), 3);
%
%       This will assign random numbers to each slice of 'mtVar' independently.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 19th November, 2010

classdef MappedTensor < handle
   properties (SetAccess = private, GetAccess = private)
      strRealFilename;        % Binary data file on disk (real part of tensor)
      strCmplxFilename;       % Binary data file on disk (complex part of tensor)
      hRealContent;           % File handle for data (real part)
      hCmplxContent;          % File handle for data (complex part)
      bTemporary;             % A flag which records whether a temporary file was created by MappedTensor
      strClass = 'double';    % The class of this mapped tensor
      strStorageClass;        % The storage class of this tensor on disk
      nClassSize;             % The size of a single scalar element of the storage class, in bytes
      vnDimensionOrder;       % A vector containing the virtual dimension order used for referencing the tensor
      nNumElements;           % The number of total elements in the tensor, for convenience
      vnOriginalSize;         % A vector recording the original size of the tensor
      bMustCast;              % A boolean indicating that the data should be cast on reading and writing
      bIsComplex = false;     % A boolean indicating the the data has a complex part
      fComplexFactor = 1;     % A factor multiplied by the complex part of the tensor (used for scalar multiplication and negation)
      fRealFactor = 1;        % A factor multiplied by the real part of the tensor (used for scalar multiplication and negation)
      nHeaderBytes = 0;       % The number of bytes to skip at the beginning of the file
   end
   
   methods
      %% MappedTensor - CONSTRUCTOR
      function [mtVar] = MappedTensor(varargin)
         % - Filter arguments for properties
         vbKeepArg = true(numel(varargin), 1);
         nArg = 2;
         while (nArg < numel(varargin))
            if (ischar(varargin{nArg}))
               switch(lower(varargin{nArg}))
                  case {'class'}
                     % - A non-default class was specified
                     mtVar.strClass = varargin{nArg+1};
                     vbKeepArg(nArg:nArg+1) = false;
                     nArg = nArg + 1;
                     
                  case {'headerbytes'}
                     % - A number of header bytes was specified
                     mtVar.nHeaderBytes = varargin{nArg+1};
                     vbKeepArg(nArg:nArg+1) = false;
                     nArg = nArg + 1;
                     
                  otherwise
                     % - No other properties are supported
                     error('MappedTensor:InvalidProperty', ...
                        '*** MappedTensor: ''%s'' is not a valid property.', varargin{nArg});
               end
            end
            
            % - Check the next argument
            nArg = nArg + 1;
         end
         
         % - Filter out unneeded arguments
         varargin = varargin(vbKeepArg);
         
         % - Get class information
         [mtVar.nClassSize, mtVar.strStorageClass] = ClassSize(mtVar.strClass);
         
         % - Do we need to cast data between these two classes?
         if (~isequal(mtVar.strStorageClass, mtVar.strClass))
            mtVar.bMustCast = true;
         else
            mtVar.bMustCast = false;
         end

         % - Should we map a file on disk, or create a temporary file?
         if (ischar(varargin{1}))
            % - Open an existing file
            vnTensorSize = [varargin{2:end}];
            mtVar.strRealFilename = varargin{1};
            mtVar.bTemporary = false;
            
         else
            % - Create a temporary file
            mtVar.bTemporary = true;
            vnTensorSize = [varargin{:}];
            
            % - Make enough space for a tensor
            mtVar.strRealFilename = CreateTempFile(prod(vnTensorSize) * mtVar.nClassSize + mtVar.nHeaderBytes);
         end
         
         % - Open the file
         mtVar.hRealContent = fopen(mtVar.strRealFilename, 'r+');
         
         % - Initialise dimension order
         mtVar.vnDimensionOrder = 1:numel(vnTensorSize);
         
         % - Record number of total elements
         mtVar.nNumElements = prod(vnTensorSize);
         
         % - Record the original tensor size
         mtVar.vnOriginalSize = vnTensorSize;
      end
      
      % delete - DESTRUCTOR
      function delete(mtVar)
         % - Delete the file, if a temporary file was created for this variable
         try
            % - Close the file handles
            fclose(mtVar.hRealContent);
            
            if (mtVar.bIsComplex)
               fclose(mtVar.hCmplxContent);
            end

            if (mtVar.bTemporary)
               % - Really delete the temporary file, don't just put it in the trash
               strState = recycle('off');
               delete(mtVar.strRealFilename);
               recycle(strState);
            end
            
            % - Delete the complex storage tensor, if it exists
            if (mtVar.bIsComplex)
               strState = recycle('off');
               delete(mtVar.strCmplxFilename);
               recycle(strState);
            end
            
         catch mtErr
            % - Die gracefully if we couldn't delete the temporary file
            warning('MappedTensor:Destructor', ...
               '--- MappedTensor/delete: Could not delete temporary file.\n       Error: %s', mtErr.message);
         end
      end
      
      %% subsref - METHOD Overloaded subsref
      function [varargout] = subsref(mtVar, subs)
         % - More than one return argument means cell or dot referencing was
         % used
         if (nargout > 1)
            error('MappedTensor:InvalidReferencing', ...
               '*** MappedTensor: ''{}'' and ''.'' referencing methods are not supported by MappedTensor objects.');
         end
         
         % - Check reference type
         switch (subs(1).type)
            case {'.', '{}'}
               % - Unsupported referencing
               error('MappedTensor:InvalidReferencing', ...
               '*** MappedTensor: ''{}'' and ''.'' referencing methods are not supported by MappedTensor objects.');
               
            case {'()'}
               % - Call the internal subsref function
               [varargout{1}] = my_subsref(mtVar, subs);

            otherwise
               % - Unknown referencing type
               error('MappedTensor:UnknownReferenceType', ...
                  '*** MappedTensor: An unknown referencing method was used.');
         end
      end
      
      % my_subsref - Standard array referencing
      function [tfData] = my_subsref(mtVar, subs)
         % - Re-order reference indices
         nNumDims = numel(subs.subs);
         nNumTotalDims = numel(mtVar.vnDimensionOrder);
         
         % - Handle different numbers of referencing dimensions
         if (nNumDims == 1)
            % - Translate from linear refs to indices
            nNumDims = numel(mtVar.vnDimensionOrder);
            
            % - Translate colon indexing
            if (isequal(subs.subs{1}, ':'))
               subs.subs{1} = (1:numel(mtVar))';
            end
            
            % - Get equivalent subscripted indexes
            vnTensorSize = size(mtVar);
            [cIndices{1:nNumDims}] = ind2sub(vnTensorSize, subs.subs{1});
            
            % - Permute indices and convert back to linear indexing
            subs.subs{1} = sub2ind(mtVar.vnOriginalSize, cIndices{mtVar.vnDimensionOrder});
            
         elseif (nNumDims < nNumTotalDims)
            % - Assume trailing dimensions are ':'
            subs.subs(nNumDims+1:numel(mtVar.vnDimensionOrder)) = {':'};
            
            % - Inverse permute index order
            vnInvOrder(mtVar.vnDimensionOrder(1:nNumTotalDims)) = 1:nNumTotalDims;
            subs.subs = subs.subs(vnInvOrder);
            
         elseif (nNumDims == nNumTotalDims)
            % - Simply permute and access tensor
            
            % - Permute index order
            vnInvOrder(mtVar.vnDimensionOrder(1:nNumTotalDims)) = 1:nNumTotalDims;
            subs.subs = subs.subs(vnInvOrder);
            
         else % (nNumDims > nNumTotalDims)
            % - Check for non-colon references
            vbNonColon = cellfun(@(c)(~isequal(c, ':')), subs.subs);
            
            % - Check only trailing dimensions
            vbNonColon(1:nNumTotalDims) = false;
            
            % - Check trailing dimensions for non-'1' indices
            if (any(cellfun(@(c)(~isequal(c, 1)), subs.subs(vbNonColon))))
               % - This is an error
               error('MappedTensor:badsubscript', ...
                  '*** MappedTensor: Index exceeds matrix dimensions.');
            end
            
            % - Permute index order
            vnInvOrder(mtVar.vnDimensionOrder(1:nNumTotalDims)) = 1:nNumTotalDims;
            subs.subs = subs.subs(vnInvOrder);            
         end
         
         % - Reference the tensor data element
         if (mtVar.bIsComplex)
            % - Get the real and complex parts
            tfData = complex(mtVar.fRealFactor .* mt_read_data(mtVar.hRealContent, subs, mtVar.vnOriginalSize, mtVar.strClass, mtVar.nHeaderBytes), ...
                             mtVar.fComplexFactor .* mt_read_data(mtVar.hCmplxContent, subs, mtVar.vnOriginalSize, mtVar.strClass, mtVar.nHeaderBytes));
         else
            % - Just return the real part
            tfData = mtVar.fRealFactor .* mt_read_data(mtVar.hRealContent, subs, mtVar.vnOriginalSize, mtVar.strClass, mtVar.nHeaderBytes);
         end
         
         % - Permute dimensions
         tfData = permute(tfData, mtVar.vnDimensionOrder);
         
         % - Reshape return data to concatenate trailing dimensions (just as
         % matlab does)
         if (nNumDims < nNumTotalDims)
            cnSize = num2cell(size(tfData));
            tfData = reshape(tfData, cnSize{1:nNumDims-1}, []);
         end
         
         % - Cast data, if required
         if (mtVar.bMustCast)
            tfData = cast(tfData, mtVar.strClass);
         end
      end
      
      %% subsasgn - METHOD Overloaded subsasgn
      function [mtVar] = subsasgn(mtVar, subs, tfData)
         % - Test real/complex nature of input and current tensor
         if (~isreal(tfData))
            % - The input data is complex
            
            if (~mtVar.bIsComplex)
               % - The tensor is currently real only
               % - Test to see if we can store complex values in the desired
               % representation
               switch (mtVar.strClass)
                  case {'char', 'logical'}
                     error('MappedTensor:NoConversionComplexToClass', ...
                        '*** MappedTensor: Cannot assign complex values to a tensor of class %s.', mtVar.strClass);
               end
               
               % - Create temporary storage for the complex part of the tensor
               if (~mtVar.bTemporary)
                  warning('MappedTensor:NoPermanentComplexStorage', ...
                     '--- MappedTensor: Warning: The complex part of a tensor is always stored temporarily.');
               end
               
               % - Make enough space for a tensor
               mtVar.strCmplxFilename = CreateTempFile(mtVar.nNumElements * mtVar.nClassSize + mtVar.nHeaderBytes);
            
               % - Open the file
               mtVar.hCmplxContent = fopen(mtVar.strCmplxFilename, 'r+');
               
               % - Record that the tensor has a complex part
               mtVar.bIsComplex = true;
            end
         end
         
         % - Cast data, if required
         if (mtVar.bMustCast)
            tfData = cast(tfData, mtVar.strStorageClass);
         end
         
         % - Permute input data
         tfData = ipermute(tfData, mtVar.vnDimensionOrder);
         
         if (~isreal(tfData))
            % - Assign to both real and complex parts
            mt_write_data(mtVar.hRealContent, subs, mtVar.vnOriginalSize, mtVar.strClass, mtVar.nHeaderBytes, real(tfData) ./ mtVar.fRealFactor);
            mt_write_data(mtVar.hCmplxContent, subs, mtVar.vnOriginalSize, mtVar.strClass, mtVar.nHeaderBytes, imag(tfData) ./ mtVar.fComplexFactor);

         else
            % - Assign only real part
            mt_write_data(mtVar.hRealContent, subs, mtVar.vnOriginalSize, mtVar.strClass, mtVar.nHeaderBytes, tfData ./ mtVar.fRealFactor);
         end
      end
      
      %% Overloaded methods (size, numel, permute, ipermute, ctranspose, transpose, isreal)
      % size - METHOD Overloaded size function
      function [varargout] = size(mtVar, vnDimensions)
         % - Return the size of the tensor data element, permuted
         vnSize = mtVar.vnOriginalSize(mtVar.vnDimensionOrder);
         
         % - Return specific dimension(s)
         if (exist('vnDimensions', 'var'))
            if (~isnumeric(vnDimensions))
               error('MappedTensor:dimensionMustBePositiveInteger', ...
                  '*** MappedTensor: Dimensions argument must be a positive integer within indexing range.');
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
            
         else
            % - Output all size elements
            varargout(1:numel(vnSize)) = num2cell(vnSize);

            % - Deal out trailing dimensions as '1'
            varargout(numel(vnSize)+1:nNumArgout) = {1};
         end
      end
      
      % numel - METHOD Overloaded numel function
      function [nNumElem] = numel(mtVar, varargin)
         % - If varargin contains anything, a cell reference "{}" was attempted
         if (~isempty(varargin))
            error('MappedTensor:cellRefFromNonCell', ...
               '*** MappedTensor: Cell contents reference from non-cell obejct.');
         end
         
         % - Return the total number of elements in the tensor
         nNumElem = mtVar.nNumElements;
      end
      
      % permute - METHOD Overloaded permute function
      function [mtVar] = permute(mtVar, vnNewOrder)
         mtVar.vnDimensionOrder(1:numel(vnNewOrder)) = mtVar.vnDimensionOrder(vnNewOrder);
      end
      
      % ipermute - METHOD Overloaded ipermute function
      function [mtVar] = ipermute(mtVar, vnOldOrder)
         vnNewOrder(vnOldOrder) = 1:numel(vnOldOrder);
         mtVar = permute(mtVar, vnNewOrder);
      end
      
      % ctranspose - METHOD Overloaded ctranspose function
      function [mtVar] = ctranspose(mtVar)
         % - Array-transpose real and complex parts
         mtVar = transpose(mtVar);
         
         % - Negate complex part
         if (mtVar.bIsComplex)
            mtVar.fComplexFactor = -mtVar.fComplexFactor;
         end
      end
      
      % transpose - METHOD Overloaded transpose function
      function [mtVar] = transpose(mtVar)
         mtVar = permute(mtVar, [2 1]);
      end
      
      % isreal - METHOD Overloaded isreal function
      function [bIsReal] = isreal(mtVar)
         bIsReal = ~mtVar.bIsComplex;
      end
      
      %% Overloaded methods (uminus, uplus, times, mtimes, ldivide, rdivide, mldivide, mrdivide)
      
      % uminus - METHOD Overloaded uminus operator (-mtVar)
      function [mtVar] = uminus(mtVar)
         % - Negate real part
         mtVar.fRealFactor = -mtVar.fRealFactor;
         
         % - Negate complex part, if it exists
         if (mtVar.bIsComplex)
            mtVar.fComplexFactor = -mtVar.fComplexFactor;
         end
      end
      
      % uplus - METHOD Overloaded uplus operator (+mtVar)
      function [mtVar] = uplus(mtVar)
         % - ...nothing to do?
      end
      
      % times - METHOD Overloaded times operator (A.*B)
      function [mtVar] = times(varargin)
         % - Are the inputs numeric?
         vbIsNumeric = cellfun(@isnumeric, varargin);
         
         % - Are the inputs scalar?
         vbIsScalar = cellfun(@isscalar, varargin);
         
         % - Are the inputs of class MappedTensor
         vbIsTensor = cellfun(@(o)(isa(o, 'MappedTensor')), varargin);
         
         % - Can we perform the operation?
         if (nnz(vbIsNumeric & vbIsScalar & ~vbIsTensor) ~= 1)
            error('MappedTensor:InvalidTimesOperands', ...
                  '*** MappedTensor: ''times'' (.*) is only supported for a MappedTensor object and a scalar.');
         end
         
         % - Get the scalar value
         fScalar = varargin{vbIsNumeric & vbIsScalar};
         mtVar = varargin{vbIsTensor};
         
         % - Multiply real and complex factors by input scalar
         mtVar.fRealFactor = mtVar.fRealFactor .* fScalar;
         
         if (mtVar.bIsComplex)
            mtVar.fCmplxFactor = mtVar.fCmplxFactor .* fScalar;
         end
      end
      
      % mtimes - METHOD Overloaded mtimes operator (A*B)
      function [mtVar] = mtimes(varargin)
         % - Use 'times' to perform the (scalar) multiplicaton
         mtVar = times(varargin{:});
      end
      
      % rdivide - METHOD Overloaded rdivide operator (A./B)
      function [mtVar] = rdivide(mtVar, fScalar)
         % - Check that the right operand is a scalar
         if (~isnumeric(fScalar) || ~isscalar(fScalar))
            error('MappedTensor:InvalidRDivideOperands', ...
                  '*** MappedTensor: ''rdivide'' (./) is only supported if the second operand is a scalar number.');
         end
         
         % - Use 'times' to perform the operation
         mtVar = times(mtVar, 1./fScalar);
      end
      
      % ldivide - METHOD Overloaded ldivide operator (A.\B)
      function [mtVar] = ldivide(fScalar, mtVar)
         % - Check that the left operand is a scalar
         if (~isnumeric(fScalar) || ~isscalar(fScalar))
            error('MappedTensor:InvalidLDivideOperands', ...
                  '*** MappedTensor: ''ldivide'' (.\\) is only supported if the first operand is a scalar number.');
         end
         
         % - Use 'times' to perform the operation
         mtVar = times(mtVar, 1./fScalar);
      end
      
      % mrdivide - METHOD Overloaded mrdivide operator (A/B)
      function [mtVar] = mrdivide(varargin)
         % - Use 'rdivide' to perform the operation
         mtVar = rdivide(varargin{:});
      end
      
      % mldivide - METHOD Overloaded mldivide operator (A\B)
      function [mtVar] = mldivide(varargin)
         % - Use 'ldivide' to perform the operation
         mtVar = ldivide(varargin{:});
      end
      
      %% Overloaded methods (plus, minus)
      
      % plus - METHOD Overloaded binary 'plus' operator (A+B)
      function [mtVar] = plus(varargin)
         % - Are the inputs numeric?
         vbIsNumeric = cellfun(@isnumeric, varargin);
         
         % - Are the inputs scalar?
         vbIsScalar = cellfun(@isscalar, varargin);
         
         % - Are the inputs of class MappedTensor
         vbIsTensor = cellfun(@(o)(isa(o, 'MappedTensor')), varargin);
         
         % - Can we perform the operation?
         if (nnz(vbIsNumeric & vbIsScalar & ~vbIsTensor) ~= 1)
            error('MappedTensor:InvalidPlusOperands', ...
                  '*** MappedTensor: ''plus'' (A+B) is only supported for a MappedTensor object and a scalar.');
         end
         
         % - Get the scalar value
         fScalar = varargin{vbIsNumeric & vbIsScalar};
         mtVar = varargin{vbIsTensor};
         
         % - Find a dimension to slice along
         vnTensorSize = size(mtVar);
         [nul, nSliceDim] = max(vnTensorSize);

         % - Perform the addition slice-wise
         fhAddition = @(tSlice)(tSlice + fScalar);
         SliceFunction(mtVar, fhAddition, nSliceDim);
      end
      
      % minus - METHOD Overloaded binary 'minus' operator (A-B)
      function [mtVar] = minus(varargin)
         % - Are the inputs numeric?
         vbIsNumeric = cellfun(@isnumeric, varargin);
         
         % - Are the inputs scalar?
         vbIsScalar = cellfun(@isscalar, varargin);
         
         % - Are the inputs of class MappedTensor
         vbIsTensor = cellfun(@(o)(isa(o, 'MappedTensor')), varargin);
         
         % - Can we perform the operation?
         if (nnz(vbIsNumeric & vbIsScalar & ~vbIsTensor) ~= 1)
            error('MappedTensor:InvalidMinusOperands', ...
                  '*** MappedTensor: ''minus'' (A-B) is only supported for a MappedTensor object and a scalar.');
         end
         
         % - Get the scalar value
         fScalar = varargin{vbIsNumeric & vbIsScalar};
         mtVar = varargin{vbIsTensor};
         
         % - Find a dimension to slice along
         vnTensorSize = size(mtVar);
         [nul, nSliceDim] = max(vnTensorSize);

         % - Perform the subtraction slice-wise
         if (vbIsScalar(1))
            fhSubtraction = @(tSlice)(fScalar - tSlice);
         else
            fhSubtraction = @(tSlice)(tSlice - fScalar);
         end
         SliceFunction(mtVar, fhSubtraction, nSliceDim);
      end
      
      %% disp - METHOD Overloaded disp function
      function disp(mtVar)
         strSize = strtrim(sprintf(' %d', size(mtVar)));
         
         if (mtVar.bIsComplex)
            strComplex = 'complex ';
         else
            strComplex = '';
         end
         
         disp(sprintf('  <a href="matlab:help MappedTensor">MappedTensor</a> class, containing: %s%s [%s].', strComplex, mtVar.strClass, strSize)); %#ok<DSPS>
         disp('  <a href="matlab:methods(''MappedTensor'')">Methods</a>');
         fprintf(1, '\n');
      end
      
      %% horzcat, vertcat, cat - METHOD Overloaded concatenation functions (unsupported)
      function out = horzcat(varargin) %#ok<STOUT,VANUS>
         error('MappedTensor:UnsupportedConcatenation', ...
            '*** MappedTensor: Concatenation is not supported for MappedTensor objects.');
      end
      function out = vertcat(varargin) %#ok<VANUS,STOUT>
         error('MappedTensor:UnsupportedConcatenation', ...
            '*** MappedTensor: Concatenation is not supported for MappedTensor objects.');
      end
      function out = cat(varargin) %#ok<VANUS,STOUT>
         error('MappedTensor:UnsupportedConcatenation', ...
            '*** MappedTensor: Concatenation is not supported for MappedTensor objects.');
      end      
   
      %% sum - METHOD Overloaded sum function for usage "sum(mtVar <, dim>)"
      function [tFinalSum] = sum(mtVar, varargin)
         % - Get tensor size
         vnTensorSize = size(mtVar);
         
         if (exist('varargin', 'var') && ~isempty(varargin))
            % - Check varargin for string parameters and discard
            vbIsString = cellfun(@ischar, varargin);
            varargin = varargin(~vbIsString);
            
            % - Too many arguments?
            if (numel(varargin) > 1)
               error('MappedTensor:sum:InvalidArguments', ...
                  '*** MappedTensor/sum: Too many arguments were supplied.');
            end
            
            % - Was a dimension specified?
            if (~isnumeric(varargin{1}) || numel(varargin{1}) > 1)
               error('MappedTensor:sum:InvalidArguments', ...
                  '*** MappedTensor/sum: ''dim'' must be supplied as a scalar number.');
            end
            
            % - Record dimension to sum along
            nDim = varargin{1};
            
         else
            % - By default, sum along first non-singleton dimension
            nDim = find(vnTensorSize > 1, 1, 'first');
         end
         
         % -- Sum in chunks to avoid allocating full tensor
         nElementsInChunk = 100000;
         vnSumSize = vnTensorSize;
         vnSumSize(nDim) = 1;
         vnSliceDimensions = cumprod(vnTensorSize);
         
         % - Find which dimension to split over (ignoring sum dimension)
         vnSplitDim = find(vnSliceDimensions >= nElementsInChunk, 2, 'first');
         nSplitDim = vnSplitDim(find(vnSplitDim ~= nDim, 1, 'first'));
         
         % - Compute the size of a single split
         vnSingleSplitSize = ceil(vnTensorSize ./ ceil(vnSliceDimensions ./ nElementsInChunk));
         
         % - Make vectors of split indices
         cellSplitIndices = cell(1, numel(vnTensorSize));
         for (nDimIndex = 1:numel(vnTensorSize)); %#ok<FORPF>
            vnStarts = 1:vnSingleSplitSize(nDimIndex):vnTensorSize(nDimIndex);
            vnEnds = [vnStarts(2:end)-1 vnTensorSize(nDimIndex)];
            vnNumDivisions(nDimIndex) = numel(vnStarts); %#ok<AGROW>
            
            for (nDivIndex = 1:vnNumDivisions(nDimIndex))
               cellSplitIndices{nDimIndex}{nDivIndex} = vnStarts(nDivIndex):vnEnds(nDivIndex);
            end
         end
         
         % -- Perform sum by taking dimensions in turn
         tFinalSum = zeros(vnSumSize);
         
         % - Construct referencing structures
         sSourceRef = substruct('()', ':');
         sDestRef = substruct('()', ':');
         
         vnSplitIndices = ones(1, numel(vnTensorSize));
         cellTheseSourceIndices = cell(1, numel(vnTensorSize));
         bContinue = true;
         while (bContinue)
            % - Find what the indices for the current chunk should be
            for (nDimIndex = 1:numel(vnTensorSize)) %#ok<FORPF>
               cellTheseSourceIndices{nDimIndex} = cellSplitIndices{nDimIndex}{vnSplitIndices(nDimIndex)};
            end
            cellTheseDestIndices = cellTheseSourceIndices;
            cellTheseDestIndices{nDim} = 1;
            
            % - Call subsasgn, subsref and sum to process data
            sSourceRef.subs = cellTheseSourceIndices;
            sDestRef.subs = cellTheseDestIndices;
            tFinalSum = subsasgn(tFinalSum, sDestRef, subsref(tFinalSum, sDestRef) + sum(subsref(mtVar, sSourceRef), nDim));
            
            % - Increment first non-max index
            nIncrementDim = find(vnSplitIndices <= vnNumDivisions, 1, 'first');
            
            % - Increment and roll-over indices, if required
            while (bContinue)
               % - Increment the index
               vnSplitIndices(nIncrementDim) = vnSplitIndices(nIncrementDim) + 1;
               
               if (vnSplitIndices(nIncrementDim) > vnNumDivisions(nIncrementDim))
                  % - We need to roll-over this index, and increment the next
                  vnSplitIndices(nIncrementDim) = 1;
                  nIncrementDim = nIncrementDim + 1;
                  
                  % - Did we roll-over the last index?
                  if (nIncrementDim > numel(vnNumDivisions))
                     bContinue = false;
                  end
                  
               else
                  % - We didn't need to roll over the index, so continue with
                  % the new indices
                  break;
               end
            end
            
         end
      end
      
      %% SliceFunction - METHOD Execute a function on the entire tensor, in slices
      function [mtNewVar] = SliceFunction(mtVar, fhFunction, nSliceDim, vnSliceSize, varargin)
         % - Get tensor size
         vnTensorSize = size(mtVar);
         
         % - Shall we generate a new tensor
         bNewTensor = false;
         
         % - Was the slice size explicity provided?
         if (~exist('vnSliceSize', 'var') || isempty(vnSliceSize))
            vnSliceSize = vnTensorSize;
            vnSliceSize(nSliceDim) = 1;
         
         elseif (~isequal(vnSliceSize([1:nSliceDim-1 nSliceDim+1:end]), vnTensorSize([1:nSliceDim-1 nSliceDim+1:end])))
            % - The slice size is different than the tensor size, so we have to
            % generate a new tensor
            bNewTensor = true;
            
            % - Display a warning if the output of this command is likely to be
            % lost
            if (nargout == 0)
               warning('MappedTensor:LostSliceOutput', ...
                  '--- MappedTensor: Warning: The output of a SliceFunction command is likely to be thrown away...');
            end
         end
         
         % - If an explicit return argument is requested, construct a new tensor
         if (nargout == 1)
            bNewTensor = true;
         end
         
         % - Shall we create a new return variable?
         if (bNewTensor)
            vnNewTensorSize = vnSliceSize;
            vnNewTensorSize(nSliceDim) = vnTensorSize(nSliceDim);
            
            mtNewVar = MappedTensor(vnNewTensorSize, 'Class', mtVar.strClass);
            
         else
            % - Store the result back in the original tensor, taking advantage
            % of the handle property of a MappedTensor
            mtNewVar = mtVar;
            
            % - Are we attempting to re-size the tensor?
            if (~isequal(vnSliceSize([1:nSliceDim-1 nSliceDim+1:end]), vnTensorSize([1:nSliceDim-1 nSliceDim+1:end])))
               error('MappedTensor:IncorrectSliceDimensions', ...
                  '*** MappedTensor/SliceFunction: A tensor cannot resized during a slice operation.\n       Assign the output to a new tensor.');
            end
         end
         
         % - Create a referencing structure
         sSubs = substruct('()', repmat({':'}, 1, numel(vnTensorSize)));
         
         % - Slice up along specified dimension
         for (nIndex = 1:vnTensorSize(nSliceDim))
            % - Make referencing structure
            sSubs.subs{nSliceDim} = nIndex;
            
            % - Handle a "slice assign" function with no input arguments efficiently
            if (nargin(fhFunction) == 0)
               subsasgn(mtNewVar, sSubs, fhFunction());
            else
               % - Execute function and store results
               subsasgn(mtNewVar, sSubs, fhFunction(subsref(mtVar, sSubs), varargin{:}));
            end
         end
      end
      
      %% saveobj - METHOD Overloaded save mechanism
      function [sVar] = saveobj(mtVar)
         % - Generate a structure containing the pertinent properties
         sVar.strRealFilename = mtVar.strRealFilename;
         sVar.bTemporary = mtVar.bTemporary;
         sVar.strClass = mtVar.strClass;
         sVar.vnDimensionOrder = mtVar.vnDimensionOrder;
         sVar.vnOriginalSize = mtVar.vnOriginalSize;
         
         % - Send a warning about poorly-supported loading
         warning('MappedTensor:UnsupportedObjectStorage', ...
            '--- MappedTensor: Warning: Saving and loaded MappedTensor objects does not preserve object data!');
      end
      
   end
      
   methods (Static)
      %% loadobj - METHOD Overloaded load mechanism
      function [mtVar] = loadobj(sSavedVar)
         % - Try to create a new MappedTensor, with the saved parameters
         if (sSavedVar.bTemporary)
            % - Create a transient mapped tensor
            mtVar = MappedTensor(sSavedVar.vnOriginalSize, 'Class', sSavedVar.strClass);
            
         else
            % - Map an existing file on disk
            mtVar = MappedTensor(sSavedVar.strRealFilename, sSavedVar.vnOriginalSize, 'Class', sSavedVar.strClass);
         end
         
         % - Record permutation
         mtVar.vnDimensionOrder = sSavedVar.vnDimensionOrder;
         
         % - Send a warning about poorly-supported loading
         warning('MappedTensor:UnsupportedObjectStorage', ...
            '--- MappedTensor: Warning: Saving and loaded MappedTensor objects does not preserve object data!');
      end
   end
end

%% --- Helper functions ---

function strFilename = CreateTempFile(nNumEntries)
   % - Get the name of a temporary file
   strFilename = tempname;
   
   % - Create the file
   hFile = fopen(strFilename, 'w+');
   
   % - Allocate enough space
   fwrite(hFile, 0, 'uint8', nNumEntries-1);
   fclose(hFile);
end

function [nBytes, strStorageClass] = ClassSize(strClass)
   % - By default, the data storage class is identical to the definition class
   strStorageClass = strClass;
   
   % - Parse class argument
   switch(lower(strClass))
      case {'char'}
         nBytes = 2;
         strStorageClass = 'uint16';
         
      case {'int8', 'uint8'}
         nBytes = 1;
         
      case {'logical'}
         nBytes = 1;
         strStorageClass = 'uint8';
         
      case {'int16', 'uint16'}
         nBytes = 2;
         
      case {'int32', 'uint32', 'single'}
         nBytes = 4;
         
      case {'int64', 'uint64', 'double'}
         nBytes = 8;
         
      otherwise
         error('MappedTensor:InvalidClass', '*** MappedTensor/ClassSize: Invalid class specifier.');
   end
end

%% Read / write functions

function [tData] = mt_read_data(hDataFile, sSubs, vnTensorSize, strClass, nHeaderBytes)
   % - Check referencing and convert to linear indices
   [vnLinearIndices, vnDataSize] = ConvertColonsCheckLims(sSubs.subs, vnTensorSize);
   
   % - Split into readable chunks
   cvnFileChunkIndices = SplitFileChunks(vnLinearIndices);
   nNumChunks = numel(cvnFileChunkIndices);
   
   % - Allocate data
   [nClassSize, strStorageClass] = ClassSize(strClass);
   tData = zeros(vnDataSize, strStorageClass);
   
   % - Read data in chunks
   nDataPointer = 1;
   for (nChunkIndex = 1:nNumChunks)
      % - Seek file to beginning of chunk
      fseek(hDataFile, (cvnFileChunkIndices{nChunkIndex}(1)-1) * nClassSize + nHeaderBytes, 'bof');
      
      % - Read chunk into return tensor
      nChunkSize = numel(cvnFileChunkIndices{nChunkIndex});
      tData(nDataPointer:nDataPointer+nChunkSize-1) = fread(hDataFile, nChunkSize, [strStorageClass '=>' strClass]);
      
      % - Shift to next data chunk
      nDataPointer = nDataPointer + nChunkSize;
   end
end

function mt_write_data(hDataFile, sSubs, vnTensorSize, strClass, nHeaderBytes, tData)
   % - Check referencing and convert to linear indices
   [vnLinearIndices, vnDataSize] = ConvertColonsCheckLims(sSubs.subs, vnTensorSize);
   
   % - Split into writable chunks
   cvnFileChunkIndices = SplitFileChunks(vnLinearIndices);
   nNumChunks = numel(cvnFileChunkIndices);

   % - Do we need to replicate the data?
   if (isscalar(tData) && prod(vnDataSize) > 1)
      tData = repmat(tData, prod(vnDataSize), 1);

   elseif (numel(tData) ~= prod(vnDataSize))
      % - The was a mismatch in the sizes of the left and right sides
      error('MappedTensor:index_assign_element_count_mismatch', ...
            '*** MappedTensor: In an assignment A(I) = B, the number of elements in B and I must be the same.');
   end
   
   % - Write data in chunks
   nDataPointer = 1;
   [nClassSize, strStorageClass] = ClassSize(strClass);
   for (nChunkIndex = 1:nNumChunks)
      % - Seek file to beginning of chunk
      fseek(hDataFile, (cvnFileChunkIndices{nChunkIndex}(1)-1) * nClassSize + nHeaderBytes, 'bof');
      
      % - Write chunk
      nChunkSize = numel(cvnFileChunkIndices{nChunkIndex});
      fwrite(hDataFile, tData(nDataPointer:nDataPointer+nChunkSize-1), strStorageClass);
      
      % - Shift to next data chunk
      nDataPointer = nDataPointer + nChunkSize;
   end
end


function [vnLinearIndices, vnDataSize] = ConvertColonsCheckLims(cRefs, vnLims)
   % - Handle linear indexing
   if (numel(cRefs) == 1)
      vnLims = prod(vnLims);
   end

   % - Check each dimension in turn
   for (nRefDim = numel(cRefs):-1:1) %#ok<FORPF>
      % - Convert colon references
      if (ischar(cRefs{nRefDim}) && isequal(cRefs{nRefDim}, ':'))
         cCheckedRefs{nRefDim} = 1:vnLims(nRefDim); %#ok<AGROW>
         
      elseif (islogical(cRefs{nRefDim}))
         % - Logical referencing -- convert to indexed referencing
         vnIndices = find(cRefs{nRefDim}(:));
         if (any(vnIndices > vnLims(nRefDim)))
            error('FocusStack:InvalidRef', ...
               '*** FocusStack/GetFullFileRefs: Logical referencing for dimension [%d] was out of bounds [1..%d].', ...
               nRefDim, vnLims(nRefDim));
         end
         cCheckedRefs{nRefDim} = vnIndices; %#ok<AGROW>
         
      elseif (any(cRefs{nRefDim}(:) < 1) || any(cRefs{nRefDim}(:) > vnLims(nRefDim)))
         % - Check limits
         error('MappedTensor:InvalidRef', ...
            '*** MappedTensor: Reference dimension [%d] was out of bounds [1..%d].', ...
            nRefDim, vnLims(nRefDim));
         
      else
         % - This dimension was ok
         cCheckedRefs{nRefDim} = cRefs{nRefDim}(:); %#ok<AGROW>
      end
   end
   
   % - Convert to linear indexing
   if (numel(cRefs) == 1)
      vnLinearIndices = cCheckedRefs{1};
   else
      [cFullRefs{1:numel(cRefs)}] = ndgrid(cCheckedRefs{:});
      vnLinearIndices = sub2ind(vnLims, cFullRefs{:});
   end
   
   % - Work out data size
   vnDataSize = cellfun(@numel, cCheckedRefs);
   
   if (numel(vnDataSize) == 1)
      vnDataSize(2) = 1;
   end
end

function [cvnFileChunkIndices] = SplitFileChunks(vnLinearIndices)
   % - Find breaks
   vnBreaks = [0 find(diff(reshape(vnLinearIndices, 1, [])) ~= 1) numel(vnLinearIndices)];
  
   % - Split indices into chunks
   for (nBreak = numel(vnBreaks)-1:-1:1)
      cvnFileChunkIndices{nBreak} = vnLinearIndices(vnBreaks(nBreak)+1:vnBreaks(nBreak+1)); %#ok<AGROW>
   end
end
% --- END of MappedTensor CLASS ---
