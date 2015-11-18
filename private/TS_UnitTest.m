% TS_UnitTest - FUNCITON Unit test for TIFFStack
%
% Usage: TS_UnitTest
%
% This function performs a series of unit tests for referencing a TIFFStack
% object, in comparison with a matlab tensor. If all tests are passed, a
% simple message reporting success will be displayed. Any error indicates a
% failure in a unit test, which will be reported in detail. Any error
% should be considered a serious bug.

% Author: Dylan Muir <dylan.muir@unibas.ch>
% Created: 17th Novemeber, 2015

function TS_UnitTest(strFilename)
   % - Make sure TIFFStack is on the path
   strTSDir = fileparts(which('TIFFStack'));
   strBaseDir = fileparts(strTSDir);
   addpath(strBaseDir);

   % - If we are being called from within 'private', then make a temp copy and run unit tests from there
   if (isequal(mfilename, 'TS_UnitTest'))
      % - Get a temporary directory and name
      strUTTestDir = fullfile(tempdir, tempname);
      mkdir(strUTTestDir);
      [~, strUTFilename] = fileparts(tempname);
      
      % - Copy the unit test file to the temporary directory
      copyfile(fullfile(strTSDir, 'private', 'TS_UnitTest.m'), fullfile(strUTTestDir, [strUTFilename '.m']));
      addpath(strUTTestDir);
      
      % - Call the unit test function
      fhUnitTest = str2func(strUTFilename);
      fhUnitTest();
      
      % - Clear up and return
      delete(fullfile(strUTTestDir, [strUTFilename '.m']));
      w = warning('off', 'MATLAB:RMDIR:RemovedFromPath');
      rmdir(strUTTestDir);
      warning(w);
      return;
   end
   
   % - Open a test TIFF file
   if (~exist('strFilename', 'var') || isempty(strFilename))
      % - Find TIFFStack and unit test file
      strFilename = fullfile(strTSDir, 'private', 'TS_UnitTestImage.tif');
   end
   
   % - Try to make a TIFFStack object
   tsStack = TIFFStack(strFilename);
   
   % - Load test image using 'imread'
   sInfo = imfinfo(strFilename);
   nNumFrames = numel(sInfo);
   
   for (nFrame = nNumFrames:-1:1)
      tfFrame = permute(imread(strFilename, 'Index', nFrame, 'Info', sInfo), [1 2 4 3]);
      tfStack(:, :, nFrame, :) = tfFrame;
   end
   
   % - Check that test stack is big enough
   vnStackSize = size(tfStack);
   if (numel(vnStackSize) < 3)
      error('TIFFStack:Arguments', '*** TS_UnitTest: Stack provided for testing [%s] must have at least 3 dimensions.', strFilename);
   end
   
   if (any(vnStackSize(1:3) < [3 3 3]))
      error('TIFFStack:Arguments', '*** TS_UnitTest: Stack provided for testing [%s] must be at least [3x3x3].', strFilename);
   end      
   
   
   %% - Test straightforward referencing
   TSUT_TestReferencing(tsStack, tfStack, 'Raw stack');
   
   %% - Test permuted stack
   tsStack = permute(tsStack, [3 1 2]);
   tfStack = permute(tfStack, [3 1 2]);
   TSUT_TestReferencing(tsStack, tfStack, 'Simple permutation');
   tsStack = ipermute(tsStack, [3 1 2]);
   tfStack = ipermute(tfStack, [3 1 2]);   
   
   %% - Test inverted stack
   tsStack.bInvert = true;
   tfStack = 255 - tfStack;
   TSUT_TestReferencing(tsStack, tfStack, 'Inverted stack');
   tsStack.bInvert = false;
   tfStack = 255 - tfStack;
   
   %% - Test non-accelerated reading funcitons
   w = warning('off', 'TIFFStack:SlowAccess');
   tsStack = TIFFStack(strFilename, [], [], true);
   TSUT_TestReferencing(tsStack, tfStack, 'Non-accelerated TIFF reading');
   warning(w);

   %% - Test deinterleaved stack
   nNumFrames = size(tfStack, 3);
   tsStack = TIFFStack(strFilename, [], [1 1 nNumFrames]);
   tfStack = reshape(tfStack, [size(tfStack, 1) size(tfStack, 2) 1 1 nNumFrames size(tfStack, 4)]);
   TSUT_TestReferencing(tsStack, tfStack, 'Deinterleaved stack');

   %% - Success if we reach here with no errors
   disp('--- TS_UnitTest: Unit tests for ''TIFFStack'' passed.');
end


function TSUT_TestReferencing(tsStack, tfStack, strTestName)
   % - Test referencing entire stack
   TSUT_compareRef(':');
   TSUT_compareRef(':', ':');
   TSUT_compareRef(':', ':', ':');
   TSUT_compareRef(':', ':', ':', ':');
   TSUT_compareRef(':', ':', ':', ':', ':');
   TSUT_compareRef(':', ':', ':', ':', ':', ':');
   
   % - Test referencing subpixels
   TSUT_compareRef(1:2, ':');
   TSUT_compareRef(1:2, 1:2, ':');
   TSUT_compareRef(1:2, 9);
   TSUT_compareRef([1 3], 9);
   TSUT_compareRef(':', ':', 3);
   
   % - Test empty referencing
   TSUT_compareRef([]);
   TSUT_compareRef([], []);
   TSUT_compareRef([], [], []);
   TSUT_compareRef([], [], [], []);
   TSUT_compareRef([], [], [], [], []);
   TSUT_compareRef(1, 1, 1, [], 1);
   
   % - Test referencing using 'end'
   TSUT_compareEndRefs();
   
   % - Test referencing with unit indices
   TSUT_compareRef(1);
   TSUT_compareRef(1, 1);
   TSUT_compareRef(1, 1, 1);
   TSUT_compareRef(1, 1, 1, 1);
   TSUT_compareRef(1, 1, 1, 1, 1);
   TSUT_compareRef(1, 1, 1, 1, 1, 1);
   
   % - Test logical indexing
   TSUT_compareRef([true false true], [true true false], [false false true]);

   % - Test limit cases that should cause errors
   vnStackSize = size(tfStack);
   TSUT_testInvalidRef('TIFFStack:badsubscript', 0);
   TSUT_testInvalidRef('TIFFStack:badsubscript', 0.1);
   TSUT_testInvalidRef('TIFFStack:badsubscript', {1});
   TSUT_testInvalidRef('TIFFStack:badsubscript', complex(1, 0));
   TSUT_testInvalidRef('TIFFStack:badsubscript', 'a');
   TSUT_testInvalidRef('TIFFStack:badsubscript', prod(vnStackSize)+1);
   TSUT_testInvalidRef('TIFFStack:badsubscript', 1, prod(vnStackSize(2:end))+1);
   TSUT_testInvalidRef('TIFFStack:badsubscript', 1, vnStackSize(2)+1, 1);
   TSUT_testInvalidRef('TIFFStack:badsubscript', 1, 1, vnStackSize(3)+1, 1);
   TSUT_testInvalidRef('TIFFStack:badsubscript', {1});
   TSUT_testInvalidRef('TIFFStack:badsubscript', 1, 1, 1, 1, 1, 1, 1, 2);
   
   % - Function to compare results of referencing stack
   function TSUT_compareRef(varargin)
      try
         % - Reference stacks, check for equality
         s.type = '()';
         s.subs = varargin;
         tfTS = subsref(tsStack, s);
         tfTF = tfStack(varargin{:});
         assert(isequal(tfTS, tfTF), 'TIFFStack:UnitTestFailed:ResultsNotEqual', ...
            'The referencing result was not equal between TIFFStack and standard tensor, with subs %s.', TSUT_subs2str(varargin));
         
      catch mErr
         mUTErr = MException('TIFFStack:UnitTestFailed', 'Failure referencing stack with subs %s, during test [%s]', ...
            TSUT_subs2str(varargin), strTestName);
         mErr = mErr.addCause(mUTErr);
         rethrow(mErr);
      end
   end

   % - Function to test referencing using 'end'
   function TSUT_compareEndRefs
      try
         % - Reference stacks, check for equality
         ae(tsStack(end), ...
            tfStack(end));
         ae(tsStack(1:end, 1), ...
            tfStack(1:end, 1));
         ae(tsStack(1, 1:end), ...
            tfStack(1, 1:end));
         ae(tsStack(1, 1:end-1), ...
            tfStack(1, 1:end-1));
         ae(tsStack(1, 1, 1:end), ...
            tfStack(1, 1, 1:end));
         ae(tsStack(1, 1, 1, 1:end), ...
            tfStack(1, 1, 1, 1:end));
         ae(tsStack(1, 1, 1, 1, 1:end), ...
            tfStack(1, 1, 1, 1, 1:end));
         ae(tsStack(1, 1, 1, 1, 1, 1:end), ...
            tfStack(1, 1, 1, 1, 1, 1:end));
         
      catch mErr
         mUTErr = MException('TIFFStack:UnitTestFailed', 'Failure referencing stack with using ''end'' referencing, during test [%s]', ...
            strTestName);
         mErr = mErr.addCause(mUTErr);
         rethrow(mErr);
      end
      
      function ae(a, b)
         assert(isequal(a, b), 'TIFFStack:UnitTestFailed:ResultsNotEqual', ...
            'The referencing result was not equal between TIFFStack and standard tensor, using ''end'' referencing');
      end
   end

   % - Function to test errors with invalid references
   function TSUT_testInvalidRef(strErrorID, varargin)
      try
         % - Reference stacks, check for equality
         s.type = '()';
         s.subs = varargin;
         tfTS = subsref(tsStack, s); %#ok<NASGU>
      
      catch mErr
         % - Check whether the correct error was reported
         if ~isequal(mErr.identifier, strErrorID)
            mUTErr = MException('TIFFStack:UnitTestFailed', 'Incorrect error raised during error test with subs %s, during test [%s]. Desired error [%s], raised error [%s].', ...
               TSUT_subs2str(varargin), strTestName, strErrorID, mErr.identifier);
            mErr = mErr.addCause(mUTErr);
            rethrow(mErr);
         end
         return;
      end
      
      % - We should never get here, so raise an error
      error('TIFFStack:UnitTestFailed:ErrorNotThrown', ...
            'An error should have occurred but did not, with subs %s, during test [%s]. Desired error [%s].', ...
            TSUT_subs2str(varargin), strTestName, strErrorID);
   end
end

function strSubs = TSUT_subs2str(varargin)
   cSubsStr = cellfun(@num2str, varargin{:}, 'UniformOutput', false);
   strSubs = ['(' ...
              sprintf('[%s], ', cSubsStr{:}) ...
              sprintf('\b\b') ...
              ')'];
end