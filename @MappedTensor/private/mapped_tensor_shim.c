/* mapped_tensor_shim.c - MEX FUNCTION Fast file access shim for MappedTensor class
 *
 * Usage:
 *
 * Author: Dylan Muir <muir@hifo.uzh.ch>
 * Created: 6th July, 2012
 */

// -- Includes

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <errno.h>


// -- MEX definitions

#ifdef   DEBUG
   #define  dprintf(...)  printf(__VA_ARGS__)
#else
   #define  dprintf(...)
#endif

#ifdef   _WIN32
   #define  FILE_FORCE_BINARY "b"
   #define  bool  uint8_t
   #define  true  1
   #define  false 0        
   #define  UINT64_C c ## i64
#else
   #define  FILE_FORCE_BINARY
#endif

#ifdef MATLAB_MEX_FILE
   #include "mex.h"
   #include "io64.h"
   #define	malloc         mxMalloc
   #define	free           mxFree
   #define  printf         mexPrintf
        
   #define  errprintf(ID, ERRTEXT, ...)  do {mexPrintf(__VA_ARGS__); mexErrMsgIdAndTxt(ID, "*** mapped_tensor_shim: " ERRTEXT);} while(0)
#else
   #error   "Compilation is only supported as a MEX file."
#endif

        
// -- Case-insensitive string comparison

#ifndef  strcasecmp
	#define  strncasecmp  cmpstrni
 
   // cmpstrni - FUNCTION Case-insensitive string comparison
   int cmpstrni(const char *s1, const char *s2, size_t length) {
      int   index;
   
      for (index = 0; *s1 && *s2 && (toupper((unsigned char) *s1) == toupper((unsigned char) *s2)) && (index <= length); s1++, s2++, index++)
         ;
   
      return *s1 - *s2;
   }
#endif
        

// -- Command definitions
        
#define  CMD_OPEN_FILE     "open"
#define  CMD_CLOSE_FILE    "close"
#define  CMD_READ_CHUNKS   "read_chunks"
#define  CMD_WRITE_CHUNKS  "write_chunks"
        

// -- Helper function definitions
        
void DisplayHelp();
char *ReadStringArg(const mxArray *pmxaArg);
mxClassID ClassIDForClassName(const char *strClassName);


// -- Command function definitions

void CmdOpenFile(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CmdCloseFile(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CmdReadChunks(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void CmdWriteChunks(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// -- Endian-management functions

inline void betoh16(uint16_t *vnData, uint64_t uNumEntries);
inline void betoh32(uint32_t *vnData, uint64_t uNumEntries);
inline void betoh64(uint64_t *vnData, uint64_t uNumEntries);
inline void letoh16(uint16_t *vnData, uint64_t uNumEntries);
inline void letoh32(uint32_t *vnData, uint64_t uNumEntries);
inline void letoh64(uint64_t *vnData, uint64_t uNumEntries);
inline void htobe16(uint16_t *vnData, uint64_t uNumEntries);
inline void htobe32(uint32_t *vnData, uint64_t uNumEntries);
inline void htobe64(uint64_t *vnData, uint64_t uNumEntries);
inline void htole16(uint16_t *vnData, uint64_t uNumEntries);
inline void htole32(uint32_t *vnData, uint64_t uNumEntries);
inline void htole64(uint64_t *vnData, uint64_t uNumEntries);

void *GetEndianSwapper(int nDataElemSize, bool bBigEndian);

// - Endian test
const int16_t i = 1;
#define is_bigendian() ( (*(int8_t *)&i) == 0 )


// -- MEX file entry function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // -- Local variables
   mwSize   nLength;
   char     *strCommand;
   
   
   // -- Check arguments
   if (nrhs < 1) {
      DisplayHelp();
      mexErrMsgIdAndTxt("MappedTensor:mapped_tensor_shim:usage",
                     "Error: Incorrect usage.");
   }
   
   
   // -- Parse command string

   // - Read command string
   strCommand = ReadStringArg(prhs[0]);

   dprintf("mts: Called with command [%s]\n", strCommand);
   
   nLength = mxGetN(prhs[0]) * mxGetM(prhs[0]) + 1;
   if (strncasecmp(strCommand, CMD_OPEN_FILE, nLength) == 0) {
      CmdOpenFile(nlhs, plhs, nrhs, prhs);
      
   } else if (strncasecmp(strCommand, CMD_CLOSE_FILE, nLength) == 0) {
      CmdCloseFile(nlhs, plhs, nrhs, prhs);
      
   } else if (strncasecmp(strCommand, CMD_READ_CHUNKS, nLength) == 0) {
      CmdReadChunks(nlhs, plhs, nrhs, prhs);
      
   } else if (strncasecmp(strCommand, CMD_WRITE_CHUNKS, nLength) == 0) {
      CmdWriteChunks(nlhs, plhs, nrhs, prhs);
      
   } else {
      DisplayHelp();
      errprintf("MappedTensor:mapped_tensor_shim:InvalidCommand",
                "Invalid command.",
                "Invalid command supplied: [%s].\n", strCommand);
   }
   
   // -- Free command string
   free(strCommand);
}
 

// -- Command functions

// Usage: [hFileHandle, strDefaultMachineFormat] = mapped_tensor_shim('open', strFilename);
void CmdOpenFile(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // -- Local variables
   char     *strFilename;
   FILE     *hFile;
   uint8_t  *uReturnFileHandle;
   
   // -- Check arguments
   
   if (nrhs < 2) {
      errprintf("MappedTensor:mapped_tensor_shim:InvalidArgument",
                "Not enough arguments for command 'open'.", "");
   }
   
   if (!mxIsChar(prhs[1])) {
      errprintf("MappedTensor:mapped_tensor_shim:InvalidArgument",
                "'strFilename' must be a character array.", "");
   }
   
   // - Get filename argument
   strFilename = ReadStringArg(prhs[1]);
   
   dprintf("mts/cof: Opening file [%s]\n", strFilename);
   
   // - Try to open file for reading
   if ((hFile = fopen(strFilename, "r+" FILE_FORCE_BINARY)) == NULL) {
      errprintf("MappedTensor:mapped_tensor_shim:CouldNotOpenFile",
                "Could not open the requested file.",
                "Could not open file [%s].  Reason: [%s]\n",
                strFilename, strerror(errno));
   }
   
   
   // -- Create the return file handle
   
   plhs[0] = mxCreateNumericMatrix(1, sizeof(FILE *), mxUINT8_CLASS, 0);
   uReturnFileHandle = (uint8_t *) mxGetData(plhs[0]);
   memcpy(uReturnFileHandle, &hFile, sizeof(FILE *));
   
   dprintf("mts/cof: File handle [%p], [%d] bytes\n", hFile, sizeof(FILE *));
   
	
	// -- Create the return default machine format string
	
	if (is_bigendian()) {
      const char  *strMachineFormat = "ieee-be";
		plhs[1] = mxCreateCharMatrixFromStrings(1, &strMachineFormat);
		dprintf("mts/cof: System format is [ieee-be]\n");
	} else {
      const char  *strMachineFormat = "ieee-le";
		plhs[1] = mxCreateCharMatrixFromStrings(1, &strMachineFormat);
		dprintf("mts/cof: System format is [ieee-le]\n");
	}
	
	
   // - Free allocated filename string
   free(strFilename);
}

// Usage: mapped_tensor_shim('close', hFileHandle);
void CmdCloseFile(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // -- Local variables
   FILE     *hFile;
   
   // -- Check arguments
   
   dprintf("mts/ccf: Closing file\n");

   if (nrhs < 2) {
      errprintf("MappedTensor:mapped_tensor_shim:InvalidArgument",
                "Not enough arguments for command 'close'.", "");
   }
   
   if ((mxGetN(prhs[1]) * mxGetN(prhs[1])) > 0) {			
		// - Read the file handle
		memcpy(&hFile, mxGetData(prhs[1]), sizeof(FILE *));
		
		dprintf("mts/ccf: File handle [%p]\n", hFile);   
		
		// - Close the file
		if (fclose(hFile) == EOF) {
			errprintf("MappedTensor:mapped_tensor_shim:CouldNotCloseFile",
						 "Could not close the requested file.",
						 "Could not close file.  Reason: [%s].\n",
						 strerror(errno));
		}
	}
}

// Usage: [tfData] = mapped_tensor_shim('read_chunks', hFileHandle, mnFileChunkIndices, vnUniqueIndices, vnReverseSort, vnDataSize, strClass, nHeaderBytes, bBigEndian);
void CmdReadChunks(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // -- Local variables
   FILE        *hFile;
   mwSize      nNumChunks, nOffset, nNumDims, *vnDataSize, nDataElemSize, nNumUniqueElems;
   double      *mfFileChunkIndices, *vfUniqueIndices, *vfReverseSort, *vfDataSize;
   char        *strClass;
   int         nChunkIndex;
   mwSize      vnSubs[2];
   mxClassID   mxcDataClass;
   uint64_t    uHeaderBytes, uElementIndex, uNumDataElems, uChunkIndex, uChunkElemIndex;
   uint8_t     *uUniqueData, *uUniqueDataPtr, *tuTargetData;
	bool			bBigEndian;
	void			(*pfEndianSwap)(void *, uint64_t);
	
   
   dprintf("mts/crc: Reading chunks\n");   
   
   // -- Check arguments
   if (nrhs < 9) {
      errprintf("MappedTensor:mapped_tensor_shim:InvalidArgument",
                "Not enough arguments for command 'read_chunks'", "");
   }
   
   // - Read the file handle
	memcpy(&hFile, mxGetData(prhs[1]), sizeof(FILE *));
   
   // - Get the chunk indices
   nNumChunks = mxGetM(prhs[2]);
   mfFileChunkIndices = (void *) mxGetPr(prhs[2]);
   
	// - Get the unique data indices
   vfUniqueIndices = mxGetPr(prhs[3]);
   nNumUniqueElems = mxGetN(prhs[3]) * mxGetM(prhs[3]);
   
   // - Get the reverse sort indices
   vfReverseSort = mxGetPr(prhs[4]);

	// - Get the output data size
   vfDataSize = mxGetPr(prhs[5]);
   nNumDims = mxGetM(prhs[5]) * mxGetN(prhs[5]);
   vnDataSize = malloc(nNumDims * sizeof(mwSize));
   uNumDataElems = 1;
   for (nOffset = 0; nOffset < nNumDims; nOffset++) {
      vnDataSize[nOffset] = (mwSize) vfDataSize[nOffset];
      uNumDataElems *= vnDataSize[nOffset];
   }
   
   // - Get the data class
   strClass = ReadStringArg(prhs[6]);
   mxcDataClass = ClassIDForClassName(strClass);
   
   // - Get the number of header bytes
   uHeaderBytes = (uint64_t) *mxGetPr(prhs[7]);
   
   // -- Allocate output variable and data cache
   
   if (mxcDataClass == mxLOGICAL_CLASS) {
      plhs[0] = mxCreateLogicalArray(nNumDims, vnDataSize);
      nDataElemSize = 1;
         
   } else {
      plhs[0] = mxCreateNumericArray(nNumDims, vnDataSize, mxcDataClass, mxREAL);
      nDataElemSize = mxGetElementSize(plhs[0]);
   }
   
   tuTargetData = (uint8_t *) mxGetPr(plhs[0]);
   if ((uUniqueData = (uint8_t *) malloc(nNumUniqueElems * nDataElemSize)) == NULL) {
		errprintf("MappedTensor:mapped_tensor_shim:Memory",
					 "Could not allocate data buffer.", "");
	}
   
	// - Get the endian-ness information and an endian-swapping function
	bBigEndian = (bool) *mxGetPr(prhs[8]);
	pfEndianSwap = GetEndianSwapper(nDataElemSize, bBigEndian);
	
	
	// -- Debug info
#ifdef   DEBUG
   dprintf("mts/crc: File handle: [%p]\n", hFile);
	
	dprintf("mts/crc: Chunks:\n");
	
	for (uChunkIndex = 0; uChunkIndex < nNumChunks; uChunkIndex++) {
		vnSubs[0] = uChunkIndex;
		
		for (uChunkElemIndex = 0, vnSubs[1] = 0; uChunkElemIndex < 3; uChunkElemIndex++, vnSubs[1]++) {
			dprintf("[%d]", (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)]);
		}
		dprintf("\n");
	}
	
	dprintf("mts/crc: Unique indices:\n");
   
	for (uChunkElemIndex = 0; uChunkElemIndex < nNumUniqueElems; uChunkElemIndex++) {
		dprintf("[%d]", (int) vfUniqueIndices[uChunkElemIndex]);
      if (uChunkElemIndex > 20) {
         dprintf("...");
         break;
      }
	}
	dprintf("\n");
	
	
	dprintf("mts/crc: Reverse sort indices:\n");
   
	for (uChunkElemIndex = 0; uChunkElemIndex < nNumUniqueElems; uChunkElemIndex++) {
		dprintf("[%d]", (int) vfReverseSort[uChunkElemIndex]);
      if (uChunkElemIndex > 20) {
         dprintf("...");
         break;
      }
	}
	dprintf("\n");
	
	dprintf("mts/crc: Data size: ");
	for(uChunkElemIndex = 0; uChunkElemIndex < nNumDims; uChunkElemIndex++) {
		dprintf("[%d]", vnDataSize[uChunkElemIndex]);
	}
	dprintf("\n");
	
	dprintf("mts/crc: Data class: [%s]\n", strClass);
	
	dprintf("mts/crc: Header bytes: [%ld]\n", uHeaderBytes);

	dprintf("mts/crc: Is big endian: [%d]\n", bBigEndian);
	
	dprintf("mts/crc: Allocated read buffer: [%d] elements, [%d] bytes per elelemt\n", nNumUniqueElems, nDataElemSize);	
#endif
   
   
   // -- Loop over chunks
   
   uUniqueDataPtr = uUniqueData;
   
   for (nChunkIndex = 0; nChunkIndex < nNumChunks; nChunkIndex++) {
      
      dprintf("mts/crc: ------- Reading chunk.  UDP is [%p] (offset %ld)\n", uUniqueDataPtr, uUniqueDataPtr - uUniqueData);
      
      // -- Get chunk info
		uint64_t	uChunkStart, uChunkSkip, uChunkSize;
		
		vnSubs[0] = nChunkIndex; vnSubs[1] = 0;
		uChunkStart = (uint64_t) (mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)] - 1) * nDataElemSize + uHeaderBytes;
		
		vnSubs[1] = 1;
		uChunkSkip  = (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)];
		
		vnSubs[1] = 2;
		uChunkSize  = (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)];
      
		
		dprintf("mts/crc: Seeking to [%ld]\n", uChunkStart);
		
      // - Seek to the beginning of this chunk
      setFilePos(hFile, (fpos_T *) &uChunkStart);

      // - Read the data for this chunk
      if (uChunkSkip == 1) {
			
         int nRead;
         
			dprintf("mts/crc: Single-read chunk, [%ld] bytes\n", uChunkSize * nDataElemSize);
			
         // - Read the data using a single read
         if ((nRead = fread((void *) uUniqueDataPtr, nDataElemSize, uChunkSize, hFile)) != uChunkSize) {
            dprintf("Read %d\n", nRead);
            if (ferror(hFile)) {
               errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                         "Error on reading file.",
                         "Could not read from file.  Reason: [%s].\n",
                         strerror(errno));
               
            } else if (feof(hFile)) {
               errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                         "Unexpected end of data file found when reading.",
                         "Single-read chunk seek [%ld] bytes; skip [0] elements; read [%ld] elements (%ld bytes)\n", 
                         uChunkStart, uChunkSize, uChunkSize * nDataElemSize);
                       
            } else {
               errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                         "Unknown error on reading file.", "");
            }
         }
      
         // - Increment data pointer
         uUniqueDataPtr += nDataElemSize * uChunkSize;
         
      } else if (uChunkSkip < 5) {
         int nRead;
         dprintf("mts/crc: Consolidated skip-read: read [%ld] bytes per element, skip [%ld] bytes.\n", nDataElemSize, nDataElemSize * (uChunkSkip-1));
         
         uint8_t  *vuConsolidatedData;
         
         if ((vuConsolidatedData = (uint8_t *) malloc(nDataElemSize * uChunkSize * uChunkSkip)) == NULL) {
            errprintf("MappedTensor:mapped_tensor_shim:Memory",
                      "Could not allocate consolidated read buffer.", "");
         }
         
         if ((nRead = fread((void *) vuConsolidatedData, 1, nDataElemSize * uChunkSkip * uChunkSize, hFile)) < nDataElemSize * uChunkSkip * uChunkSize) {
            dprintf("Read [%ld] bytes; expected [%ld] bytes; limit [%ld]\n", 
                    nRead, nDataElemSize * uChunkSkip * uChunkSize,
                    nDataElemSize * uChunkSkip * (uChunkSize-1) + nDataElemSize);
            
            // - We might get a short read if we're at the end of the file
            if (nRead < nDataElemSize * uChunkSkip * (uChunkSize-1) + nDataElemSize) {  
               if (ferror(hFile)) {
                  errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                          "Error on reading file.",
                          "Could not read from file.  Reason: [%s].\n",
                          strerror(errno));
                  
               } else if (feof(hFile)) {
               errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                         "Unexpected end of data file found when reading.",
                         "Consolidated skip-read chunk seek [%ld] bytes; skip [%ld] elements; read [%ld] elements (%ld bytes)\n", 
                         uChunkStart, uChunkSkip, uChunkSize, uChunkSize * nDataElemSize);
                  
               } else {
                  errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                          "Unknown error on reading file.", "");
               }
            }
         }
         
         // - Copy into unique data buffer
         for (uElementIndex = 0; uElementIndex < uChunkSize; uElementIndex++) {
            memcpy((void *) uUniqueDataPtr + (uElementIndex * nDataElemSize), vuConsolidatedData + (uElementIndex * nDataElemSize * uChunkSkip), nDataElemSize);
         }
         
         // - Free buffer
         free(vuConsolidatedData);
         
         // - Increment data pointer
         uUniqueDataPtr = uUniqueDataPtr + uChunkSize * nDataElemSize * (uChunkSkip-1);
         
         
      } else {
      
			dprintf("mts/crc: Single-element skip-read: read [%ld] bytes per element, skip [%ld] bytes.\n", nDataElemSize, nDataElemSize * (uChunkSkip-1));
			
         // - Read an element, then skip elements
         for (uElementIndex = 0; uElementIndex < uChunkSize; uElementIndex++) {
            // - Read an element
            if (fread((void *) uUniqueDataPtr, nDataElemSize, 1, hFile) != 1) {
               if (ferror(hFile)) {
                  errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                          "Error on reading file.",
                          "Could not read from file.  Reason: [%s].\n",
                          strerror(errno));
                  
               } else if (feof(hFile)) {
               errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                         "Unexpected end of data file found when reading.",
                         "Single skip-read chunk seek [%ld] bytes; skip [%ld] elements; read [1] elements (%ld bytes)\n", 
                         uChunkStart, uChunkSkip, nDataElemSize);
                  
               } else {
                  errprintf("MappedTensor:mapped_tensor_shim:FileReadError",
                          "Unknown error on reading file.", "");
               }
            }
            
            // - Skip the required number of elements
            if (uElementIndex < (uChunkSize-1)) {
               dprintf("mts/crc: Skipping [%ld] bytes\n", (uChunkSkip-1) * nDataElemSize);
               fseek(hFile, (uChunkSkip-1) * nDataElemSize, SEEK_CUR);
            }
            
            // - Move to the next buffer element
            uUniqueDataPtr += nDataElemSize;
         }
      }
   }
   
   // - Assign unique data to target data, in original indexing order
   for (uElementIndex = 0; uElementIndex < uNumDataElems; uElementIndex++) {
      memcpy(tuTargetData + (uElementIndex * nDataElemSize), uUniqueData + ((uint64_t)(vfReverseSort[uElementIndex] - 1) * nDataElemSize), nDataElemSize);
   }
   
	// - Perform endian swap, if necessary
	if (pfEndianSwap != NULL) {
		pfEndianSwap(tuTargetData, uNumDataElems);
	}
	
   // - Free temporary storage
   free(vnDataSize);
   free(uUniqueData);
}


/*
% mt_read_data_chunks - FUNCTION Read data without sorting or checking indices
% 'vnUniqueIndices' MUST be sorted and unique; 'vnReverseSort' must be the
% inverting indices from calling UNIQUE
function [tData] = mt_read_data_chunks(hDataFile, mnFileChunkIndices, vnUniqueIndices, vnReverseSort, vnDataSize, strClass, nHeaderBytes)
   nNumChunks = size(mnFileChunkIndices, 1);
   
   % - Allocate data
   [nClassSize, strStorageClass] = ClassSize(strClass);
   vUniqueData = zeros(numel(vnUniqueIndices), 1, strStorageClass);
   
   % - Read data in chunks
   nDataPointer = 1;
   for (nChunkIndex = 1:nNumChunks)
      % - Get chunk info
      nChunkSkip = mnFileChunkIndices(nChunkIndex, 2);
      nChunkSize = mnFileChunkIndices(nChunkIndex, 3);
      
      % - Seek file to beginning of chunk
      fseek(hDataFile, (mnFileChunkIndices(nChunkIndex, 1)-1) * nClassSize + nHeaderBytes, 'bof');
      
      % - Normal forward read
      vUniqueData(nDataPointer:nDataPointer+nChunkSize-1) = fread(hDataFile, nChunkSize, [strStorageClass '=>' strClass], (nChunkSkip-1) * nClassSize);
      
      % - Shift to next data chunk
      nDataPointer = nDataPointer + nChunkSize;
   end
   
   % - Assign data back to original indexing order
   tData = reshape(vUniqueData(vnReverseSort), vnDataSize);
end

*/

// Usage: mapped_tensor_shim('write_chunks', hFileHandle, cvnFileChunkIndices, vnUniqueDataIndices, vnDataSize, strClass, nHeaderBytes, tData, bBigEndian);
void CmdWriteChunks(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // -- Local variables
   FILE        *hFile;
   mwSize      nNumChunks, nOffset, nNumDims, *vnDataSize, nDataElemSize, nNumUniqueElems;
   double      *mfFileChunkIndices, *vfUniqueIndices, *vfDataSize;
   char        *strClass;
   mxClassID   mxcDataClass;
   mxArray     *mxaTest;
   uint64_t    uHeaderBytes, uUniqueElemIndex, uChunkElemIndex, uNumRefElems, uNumDataElems, uMaxChunkSize, uChunkIndex;
   uint8_t     *uSourceData, *tuTargetData;
   bool        bScalarData, bBigEndian;
	mwSize		vnSubs[2];
	uint64_t		uChunkStart, uChunkSkip, uChunkSize;
	void			(*pfEndianSwap)(void *, uint64_t);
   
   dprintf("mts/cwc: Writing chunks\n");   
   
   // -- Check arguments
   if (nrhs < 9) {
      errprintf("MappedTensor:mapped_tensor_shim:InvalidArgument",
                "Not enough arguments for command 'write_chunks'", "");
   }
   
   // - Read the file handle
	memcpy(&hFile, mxGetData(prhs[1]), sizeof(FILE *));
   
   // - Get the chunk indices
   nNumChunks = mxGetM(prhs[2]);
   mfFileChunkIndices = mxGetPr(prhs[2]);
   
   // - Get the unique data indices
   vfUniqueIndices = mxGetPr(prhs[3]);
   nNumUniqueElems = mxGetN(prhs[3]) * mxGetM(prhs[3]);
   
   // - Get the referenced data size
   vfDataSize = mxGetPr(prhs[4]);
   nNumDims = mxGetM(prhs[4]) * mxGetN(prhs[4]);
   vnDataSize = malloc(nNumDims * sizeof(mwSize));
   uNumRefElems = 1;
   for (nOffset = 0; nOffset < nNumDims; nOffset++) {
      vnDataSize[nOffset] = (mwSize) vfDataSize[nOffset];
      uNumRefElems *= vnDataSize[nOffset];
   }
   
   // - Get the data class
   strClass = ReadStringArg(prhs[5]);
   mxcDataClass = ClassIDForClassName(strClass);
   
   // - Get the number of header bytes
   uHeaderBytes = (uint64_t) *mxGetPr(prhs[6]);
   
   // - Get the data to be written and element size
   uSourceData = (uint8_t *) mxGetPr(prhs[7]);
   uNumDataElems = mxGetM(prhs[7]) * mxGetN(prhs[7]);
   bScalarData = uNumDataElems == 1;
   
   if (mxcDataClass == mxLOGICAL_CLASS) {
      nDataElemSize = 1;
         
   } else {
      mxaTest = mxCreateNumericArray(nNumDims, vnDataSize, mxcDataClass, mxREAL);
      nDataElemSize = mxGetElementSize(mxaTest);
      mxDestroyArray(mxaTest);
   }
   
	// - Get the endian-ness information and an endian-swapping function
	bBigEndian = (bool) *mxGetPr(prhs[8]);
	pfEndianSwap = GetEndianSwapper(nDataElemSize, bBigEndian);

	
	// -- Debug info
#ifdef   DEBUG
   dprintf("mts/cwc: File handle: [%p]\n", hFile);
	
	dprintf("mts/cwc: Chunks:\n");
	
	for (uChunkIndex = 0; uChunkIndex < nNumChunks; uChunkIndex++) {
		vnSubs[0] = uChunkIndex;
		
		for (uChunkElemIndex = 0, vnSubs[1] = 0; uChunkElemIndex < 3; uChunkElemIndex++, vnSubs[1]++) {
			dprintf("[%d]", (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)]);
		}
		dprintf("\n");
	}
	
	dprintf("mts/cwc: Unique indices:\n");
   
	for (uChunkElemIndex = 0; uChunkElemIndex < nNumUniqueElems; uChunkElemIndex++) {
		dprintf("[%d]", (int) vfUniqueIndices[uChunkElemIndex]);
      if (uChunkElemIndex > 20) {
         dprintf("...");
         break;
      }
	}
	dprintf("\n");
		
	if (bScalarData) {
		dprintf("mts/cwc: Scalar source data.\n");
	} else {
		dprintf("mts/cwc: Source data size: [%ld] elements\n", uNumDataElems);
	}
	
	dprintf("mts/cwc: Reference data size: ");
	for(uChunkElemIndex = 0; uChunkElemIndex < nNumDims; uChunkElemIndex++) {
		dprintf("[%d]", vnDataSize[uChunkElemIndex]);
      if (uChunkElemIndex > 20) {
         dprintf("...");
         break;
      }
	}
	dprintf("\n");
	
	dprintf("mts/cwc: Data class: [%s]\n", strClass);
	
	dprintf("mts/cwc: Header bytes: [%ld]\n", uHeaderBytes);
	
	dprintf("mts/cwc: Is big endian: [%d]\n", bBigEndian);
#endif
	
   // - Check data and reference sizes
   if (!bScalarData && (uNumDataElems != uNumRefElems)) {
      errprintf("MappedTensor:index_assign_element_count_mismatch",
                "In an assignment A(I) = B, the number of elements in B and I must be the same.", "");
   }
	
	// -- Get max chunk size, and allocate data buffer
	
	uMaxChunkSize = 0;
   for (uChunkIndex = 0; uChunkIndex < nNumChunks; uChunkIndex++) {
      // -- Get chunk size
		vnSubs[0] = uChunkIndex;
		vnSubs[1] = 2;
      if (mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)] > uMaxChunkSize)
			uMaxChunkSize = (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)];
	}
	
	if ((tuTargetData = (uint8_t *) malloc(uMaxChunkSize * nDataElemSize)) == NULL) {
		errprintf("MappedTensor:mapped_tensor_shim:Memory",
					 "Could not allocate data buffer.", "");
	}
	
	dprintf("mts/cwc: Allocated chunk data write buffer: [%d] elements, [%d] bytes per elelemt\n", uMaxChunkSize, nDataElemSize);	
	
	
	// -- If scalar data, fill data write buffer with replicated copies of scalar data
	
	if (bScalarData) {
		dprintf("mts/cwc: Replicating scalar data to write buffer.\n");
		for (uChunkElemIndex = 0; uChunkElemIndex < uMaxChunkSize; uChunkElemIndex++) {
			memcpy((void *) (tuTargetData + uChunkElemIndex * nDataElemSize), (void *) uSourceData, nDataElemSize);
		}
	}
	
	
   // -- Write data
	
	uUniqueElemIndex = 0;
   for (uChunkIndex = 0; uChunkIndex < nNumChunks; uChunkIndex++) {
      // -- Get chunk info		
		vnSubs[0] = uChunkIndex; vnSubs[1] = 0;
		uChunkStart = (uint64_t) (mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)] - 1) * nDataElemSize + uHeaderBytes;
		
		vnSubs[1] = 1;
		uChunkSkip  = (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)];
		
		vnSubs[1] = 2;
		uChunkSize  = (uint64_t) mfFileChunkIndices[mxCalcSingleSubscript(prhs[2], 2, vnSubs)];
      
      // - Seek to the beginning of this chunk
      setFilePos(hFile, (fpos_T *) &uChunkStart);
      
      dprintf("mts/cwc: Seeking to [%ld] bytes\n", uChunkStart);
		
		// - Gather data for this chunk
		if (!bScalarData) {
			for (uChunkElemIndex = 0; uChunkElemIndex < uChunkSize; uChunkElemIndex++, uUniqueElemIndex++) {
				memcpy((void *) (tuTargetData + uChunkElemIndex * nDataElemSize), (void *) uSourceData + (((uint64_t) (vfUniqueIndices[uUniqueElemIndex]-1)) * nDataElemSize), nDataElemSize);
			}
		}
      
		// - Perform endian swap, if necessary
		if (pfEndianSwap != NULL) {
			pfEndianSwap(tuTargetData, uChunkSize);
		}
		
      // - Write data for this chunk
      if (uChunkSkip == 1) {
         // - Write the data using a single write
			dprintf("mts/cwc: Single-write chunk [%ld] bytes\n", uChunkSize * nDataElemSize);
			
         if (fwrite((void *) tuTargetData, nDataElemSize, uChunkSize, hFile) != uChunkSize) {
            if (ferror(hFile)) {
               errprintf("MappedTensor:mapped_tensor_shim:FileWriteError",
                         "Error on writing file.",
                         "Could not write to file.  Reason: [%s].\n",
                         strerror(errno));
               
            } else {
               errprintf("MappedTensor:mapped_tensor_shim:FileWriteError",
                         "Unknown error on writing to file.", "");
            }
         }
			
      } else {
			
         // - Write an element, then skip elements
			dprintf("mts/cwc: Skip-write chunk: [%ld] bytes per element\n", nDataElemSize);
			
         for (uChunkElemIndex = 0; uChunkElemIndex < uChunkSize; uChunkElemIndex++) {
            // - Write an element
            if (fwrite((void *) (tuTargetData + (uChunkElemIndex * nDataElemSize)), nDataElemSize, 1, hFile) != 1) {
					if (ferror(hFile)) {
						errprintf("MappedTensor:mapped_tensor_shim:FileWriteError",
									 "Error on writing file.",
									 "Could not write to file.  Reason: [%s].\n",
									 strerror(errno));
						
					} else {
						errprintf("MappedTensor:mapped_tensor_shim:FileWriteError",
									 "Unknown error on writing to file.", "");
					}
            }
            
            // - Skip the required number of elements (one was already skipped due to the write above)
            fseek(hFile, (uChunkSkip-1) * nDataElemSize, SEEK_CUR);
				
				dprintf("mts/cwc: Skip-write chunk: Skipping over [%ld] bytes\n", (uChunkSkip-1) * nDataElemSize);
         }
      }
   }
	
	// - Flush stream
	if (fflush(hFile) == EOF) {
		errprintf("MappedTensor:mapped_tensor_shim:FileWriteError",
					 "Error on writing file.  Could not flush stream",
					 "Could not write to file.  Reason: [%s].\n",
					 strerror(errno));
	}
	
	dprintf("mts/cwc: Flushed stream\n");
}

/*
% mt_write_data_chunks - FUNCTION Write data without sorting or checking indices
% 'vnUniqueIndices' MUST be sorted and unique; 'vnUniqueDataIndices' must
% be the corresponding indices into the data from calling UNIQUE
function mt_write_data_chunks(hDataFile, mnFileChunkIndices, vnUniqueDataIndices, vnDataSize, strClass, nHeaderBytes, tData)
   nNumChunks = size(mnFileChunkIndices, 1);

   % - Do we need to replicate the data?
   if (isscalar(tData) && prod(vnDataSize) > 1)
      tData = repmat(tData, prod(vnDataSize), 1);

   elseif (numel(tData) ~= prod(vnDataSize))
      % - The was a mismatch in the sizes of the left and right sides
      error('MappedTensor:index_assign_element_count_mismatch', ...
            '*** MappedTensor: In an assignment A(I) = B, the number of elements in B and I must be the same.');
   end
   
   % - Take only unique data indices
   vUniqueData = tData(vnUniqueDataIndices);
   
   % - Write data in chunks
   nDataPointer = 1;
   [nClassSize, strStorageClass] = ClassSize(strClass);
   for (nChunkIndex = 1:nNumChunks)
      % - Get chunk info
      nChunkSkip = mnFileChunkIndices(nChunkIndex, 2);
      nChunkSize = mnFileChunkIndices(nChunkIndex, 3);

      % - Seek file to beginning of chunk
      fseek(hDataFile, (mnFileChunkIndices(nChunkIndex, 1)-1) * nClassSize + nHeaderBytes, 'bof');
      
      % - Normal forward write of chunk data
      fwrite(hDataFile, vUniqueData(nDataPointer:nDataPointer+nChunkSize-1), strStorageClass, (nChunkSkip-1) * nClassSize);
      
      % - Shift to next data chunk
      nDataPointer = nDataPointer + nChunkSize;
   end
end
*/




// -- Helper functions

void DisplayHelp()
{
   printf("mapped_tensor_shim - MEX FUNCTION Fast file access shim for MappedTensor class\n\n");
   
   printf("   Usage: [hFileHandle] = mapped_tensor_shim('open', strFilename);\n");
   printf("          mapped_tensor_shim('close', hFileHandle);\n");
   printf("          [tfData] = mapped_tensor_shim('read_chunks', hFileHandle, mnFileChunkIndices, vnUniqueIndices, vnReverseSort, vnDataSize, strClass, nHeaderBytes);\n");
   printf("          mapped_tensor_shim('write_chunks', hFileHandle, cvnFileChunkIndices, vnUniqueDataIndices, vnDataSize, strClass, nHeaderBytes, tData);\n");
   printf("\n");
}

char *ReadStringArg(const mxArray *pmxaArg)
{
   mwSize   nLength;
   char     *strCommand;
   
   nLength = mxGetN(pmxaArg) * mxGetM(pmxaArg) + 1;
   strCommand = (char *) malloc(nLength * sizeof(mxChar));
   
   if (mxGetString(pmxaArg, strCommand, nLength) == 1) {
      errprintf("MappedTensor:mapped_tensor_shim:InvalidArgument",
                "Could not read a string argument.", "");
   }
	
	return strCommand;
}

mxClassID ClassIDForClassName(const char *strClassName)
{
   if (strcasecmp(strClassName, "cell") == 0) {
      return mxCELL_CLASS;
      
   } else if (strcasecmp(strClassName, "char") == 0) {
      return mxCHAR_CLASS;
      
   } else if (strcasecmp(strClassName, "double") == 0) {
      return mxDOUBLE_CLASS;
      
   } else if (strcasecmp(strClassName, "function_handle") == 0) {
      return mxFUNCTION_CLASS;
      
   } else if (strcasecmp(strClassName, "int8") == 0) {
      return mxINT8_CLASS;
      
   } else if (strcasecmp(strClassName, "int16") == 0) {
      return mxINT16_CLASS;
         
   } else if (strcasecmp(strClassName, "int32") == 0) {
      return mxINT32_CLASS;
      
   } else if (strcasecmp(strClassName, "int64") == 0) {
      return mxINT64_CLASS;
      
   } else if (strcasecmp(strClassName, "logical") == 0) {
      return mxLOGICAL_CLASS;
      
   } else if (strcasecmp(strClassName, "single") == 0) {
      return mxSINGLE_CLASS;
         
   } else if (strcasecmp(strClassName, "struct") == 0) {
      return mxSTRUCT_CLASS;
      
   } else if (strcasecmp(strClassName, "uint8") == 0) {
      return mxUINT8_CLASS;
      
   } else if (strcasecmp(strClassName, "uint16") == 0) {
      return mxUINT16_CLASS;
      
   } else if (strcasecmp(strClassName, "uint32") == 0) {
      return mxUINT32_CLASS;

   } else if (strcasecmp(strClassName, "uint64") == 0) {
      return mxUINT64_CLASS;
        
   } else if (strcasecmp(strClassName, "unknown") == 0) {
      return mxUNKNOWN_CLASS;
      
   } else {
      return mxUNKNOWN_CLASS;
   }
}


// --- Endian managing functions

inline uint16_t ByteSwap16(uint16_t *nData) {
   *nData = ((*nData & 0xFF) << 8) | ((*nData & 0xFF00) >> 8);
}

inline uint32_t ByteSwap32(uint32_t *nData) {
   *nData = ((*nData & 0xFF) << 24) | ((*nData & 0xFF00) << 8) | ((*nData & 0xFF0000) >> 8) | ((*nData & 0xFF000000) >> 24);
}

inline uint64_t ByteSwap64(uint64_t *nData)
{
   *nData = 
   ((*nData &   UINT64_C(0x00000000000000FF)) << 56) |
   ((*nData &   UINT64_C(0x000000000000FF00)) << 40) |
   ((*nData &   UINT64_C(0x0000000000FF0000)) << 24) |
   ((*nData &   UINT64_C(0x00000000FF000000)) << 8)  |
   ((*nData &   UINT64_C(0x000000FF00000000)) >> 8)  | 
   ((*nData &   UINT64_C(0x0000FF0000000000)) >> 24) |
   ((*nData &   UINT64_C(0x00FF000000000000)) >> 40) |
   ((*nData &   UINT64_C(0xFF00000000000000)) >> 56);
}

void ByteSwapBuffer16(uint16_t *vnData, uint64_t uNumEntries) {
   uint64_t nIndex;
   
   for (nIndex = 0; nIndex < uNumEntries; nIndex++) {
		ByteSwap16(&vnData[nIndex]);
   }
}

void ByteSwapBuffer32(uint32_t *vnData, uint64_t uNumEntries) {
   uint64_t nIndex;
   
   for (nIndex = 0; nIndex < uNumEntries; nIndex++) {
		ByteSwap32(&vnData[nIndex]);
   }
}

void ByteSwapBuffer64(uint64_t *vnData, uint64_t uNumEntries) {
   uint64_t nIndex;
   
   for (nIndex = 0; nIndex < uNumEntries; nIndex++) {
		ByteSwap64(&vnData[nIndex]);
   }
}

inline void betoh16(uint16_t *vnData, uint64_t uNumEntries) {
	if (!is_bigendian())
		ByteSwapBuffer16(vnData, uNumEntries);
}

inline void betoh32(uint32_t *vnData, uint64_t uNumEntries) {
	if (!is_bigendian())
		ByteSwapBuffer32(vnData, uNumEntries);
}

inline void betoh64(uint64_t *vnData, uint64_t uNumEntries) {
	if (!is_bigendian())
		ByteSwapBuffer64(vnData, uNumEntries);
}

inline void letoh16(uint16_t *vnData, uint64_t uNumEntries) {
	if (is_bigendian())
		ByteSwapBuffer16(vnData, uNumEntries);
}

inline void letoh32(uint32_t *vnData, uint64_t uNumEntries) {
	if (is_bigendian())
		ByteSwapBuffer32(vnData, uNumEntries);
}

inline void letoh64(uint64_t *vnData, uint64_t uNumEntries) {
	if (is_bigendian())
		ByteSwapBuffer64(vnData, uNumEntries);
}

inline void htobe16(uint16_t *vnData, uint64_t uNumEntries) {
	if (!is_bigendian())
		ByteSwapBuffer16(vnData, uNumEntries);
}

inline void htobe32(uint32_t *vnData, uint64_t uNumEntries) {
	if (!is_bigendian())
		ByteSwapBuffer32(vnData, uNumEntries);
}

inline void htobe64(uint64_t *vnData, uint64_t uNumEntries) {
	if (!is_bigendian())
		ByteSwapBuffer64(vnData, uNumEntries);
}

inline void htole16(uint16_t *vnData, uint64_t uNumEntries) {
	if (is_bigendian())
		ByteSwapBuffer16(vnData, uNumEntries);
}

inline void htole32(uint32_t *vnData, uint64_t uNumEntries) {
	if (is_bigendian())
		ByteSwapBuffer32(vnData, uNumEntries);
}

inline void htole64(uint64_t *vnData, uint64_t uNumEntries) {
	if (is_bigendian())
		ByteSwapBuffer64(vnData, uNumEntries);
}

void *GetEndianSwapper(int nDataElemSize, bool bBigEndian) {
	switch (nDataElemSize) {
		case 1:
			return NULL;
			
		case 2:
			if (bBigEndian) {
				return &htobe16;
			} else {
				return &htole16;
			}
			
		case 4:
			if (bBigEndian) {
				return &htobe32;
			} else {
				return &htole32;
			}
			
		case 8:
			if (bBigEndian) {
				return &htobe64;
			} else {
				return &htole64;
			}

		default:
			errprintf("MappedTensor:mapped_tensor_shim:EndianDataSize",
						 "Unknown data size for endian swapping.", "");
	}
}

/* --- END of mapped_tensor_shim.c ---*/
