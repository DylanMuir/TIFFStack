/******************
 * bin_frames_time - Average frames with temporal binnning
 *
 * Usage: bin_frames_time strInputFile strOutputFile nXPixels nYPixels nChannels nBytesPerPixel nFramesPerBin bBigEndian
 */

/* Author: Dylan Muir <muir@hifo.uzh.ch>
 * Created: 2nd July, 2012
 */

/* - Includes */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
/*#include <stdbool.h>*/
#define  bool  uint8_t

#ifdef   _WIN32
   #define  FORCE_BINARY_FILE "b"
#else
   #define  FORCE_BINARY_FILE
#endif

#ifdef MATLAB_MEX_FILE
   #include "mex.h"
   #define	malloc         mxMalloc
   #define	free           mxFree
   #define  printf         mexPrintf
        
   #define  errprintf(...)  do {mexPrintf(__VA_ARGS__); mexErrMsgTxt("Error (see text above)");} while(0)
#else
   #define  errprintf(...) fprintf(stderr, __VA_ARGS__)
#endif


/* - Endian test */
const int16_t i = 1;
#define is_bigendian() ( (*(int8_t *)&i) == 0 )


uint16_t ByteSwapShort(uint16_t nData) {
   return ((nData & 0xFF) << 8) | ((nData & 0xFF00) >> 8);
}

uint32_t ByteSwapLong(uint32_t nData) {
   return ((nData & 0xFF) << 24) | ((nData & 0xFF00) << 8) | ((nData & 0xFF0000) >> 8) | ((nData & 0xFF000000) >> 24);
}

void ByteSwapBufferShort(uint16_t *vnSourceData, uint64_t *vnDestData, long unsigned int nNumEntries, bool bSwap) {
   uint16_t				nThisData;
   unsigned long int nIndex;
   
   for (nIndex = 0; nIndex < nNumEntries; nIndex++) {
      if (bSwap)
         nThisData = ByteSwapShort(vnSourceData[nIndex]);
      else
         nThisData = vnSourceData[nIndex];
      
      vnDestData[nIndex] = (unsigned long int) nThisData;
   }
}

void ByteSwapBufferLong(uint32_t *vnSourceData, uint64_t *vnDestData, long unsigned int nNumEntries, bool bSwap) {
   uint32_t				nThisData;
   unsigned long int nIndex;
   
   for (nIndex = 0; nIndex < nNumEntries; nIndex++) {
      if (bSwap)
         nThisData = ByteSwapLong(vnSourceData[nIndex]);
      else
         nThisData = vnSourceData[nIndex];
      
      vnDestData[nIndex] = (unsigned long int) nThisData;
   }
}

void BufferSwapByteShort(uint64_t *vnSourceData, uint16_t *vnDestData, long unsigned int nNumEntries, bool bSwap) {
   uint16_t				nThisData;
   unsigned long int nIndex;
   
   for (nIndex = 0; nIndex < nNumEntries; nIndex++) {
      nThisData = (uint16_t) vnSourceData[nIndex];
      
      if (bSwap)
         vnDestData[nIndex] = ByteSwapShort(nThisData);
      else
         vnDestData[nIndex] = nThisData;
   }
}

void BufferSwapByteLong(uint64_t *vnSourceData, uint32_t *vnDestData, long unsigned int nNumEntries, bool bSwap) {
   uint32_t				nThisData;
   unsigned long int	nIndex;
   
   for (nIndex = 0; nIndex < nNumEntries; nIndex++) {
      nThisData = (uint32_t) vnSourceData[nIndex];
      
      if (bSwap)
         vnDestData[nIndex] = ByteSwapLong(nThisData);
      else
         vnDestData[nIndex] = nThisData;
   }
}

void betohs(uint16_t *vnSourceData, uint64_t *vnDestData, long unsigned int nNumEntries) {
   ByteSwapBufferShort(vnSourceData, vnDestData, nNumEntries, !is_bigendian());
}

void betohl(uint32_t *vnSourceData, uint64_t *vnDestData, long unsigned int nNumEntries) {
   ByteSwapBufferLong(vnSourceData, vnDestData, nNumEntries, !is_bigendian());
}

void letohs(uint16_t *vnSourceData, uint64_t *vnDestData, long unsigned int nNumEntries) {
   ByteSwapBufferShort(vnSourceData, vnDestData, nNumEntries, is_bigendian());
}

void letohl(uint32_t *vnSourceData, uint64_t *vnDestData, long unsigned int nNumEntries) {
   ByteSwapBufferLong(vnSourceData, vnDestData, nNumEntries, is_bigendian());
}

void htobes(uint64_t *vnSourceData, uint16_t *vnDestData, long unsigned int nNumEntries) {
   BufferSwapByteShort(vnSourceData, vnDestData, nNumEntries, !is_bigendian());
}

void htobel(uint64_t *vnSourceData, uint32_t *vnDestData, long unsigned int nNumEntries) {
   BufferSwapByteLong(vnSourceData, vnDestData, nNumEntries, !is_bigendian());
}

void htoles(uint64_t *vnSourceData, uint16_t *vnDestData, long unsigned int nNumEntries) {
   BufferSwapByteShort(vnSourceData, vnDestData, nNumEntries, is_bigendian());
}

void htolel(uint64_t *vnSourceData, uint32_t *vnDestData, long unsigned int nNumEntries) {
   BufferSwapByteLong(vnSourceData, vnDestData, nNumEntries, is_bigendian());
}


int ReadFrame(FILE *hFile, uint64_t *vnFrameBuffer, uint8_t *vnFrameBytesBuf, unsigned long int nPixelsPerFrame, long unsigned int nBytesPerPixel, bool bBigEndian) {
   int nRead, nIndex;
   
   if (feof(hFile))
      return -1;
   
   /* - Read a frame in bytes */
   nRead = fread((void *)vnFrameBytesBuf, 1, nPixelsPerFrame * nBytesPerPixel, hFile);
   
   if (nRead != (nPixelsPerFrame * nBytesPerPixel)) {
      if (ferror(hFile)) {
         errprintf("*** bin_frames_time/ReadFrame: Could not read from file.\n   Reason: %s\n", strerror(errno));
      }

      return -1;
      
   } else {
      /* Convert endian-ness */
      switch (nBytesPerPixel) {
         case 1:
            for (nIndex = 0; nIndex < nPixelsPerFrame; nIndex++)
               vnFrameBuffer[nIndex] = (unsigned long) vnFrameBytesBuf[nIndex];
            
            break;
            
         case 2:
            if (bBigEndian) {
               betohs((uint16_t *) vnFrameBytesBuf, vnFrameBuffer, nPixelsPerFrame);
            } else {
               letohs((uint16_t *) vnFrameBytesBuf, vnFrameBuffer, nPixelsPerFrame);
            }
            break;
            
         case 4:
            if (bBigEndian) {
               betohl((uint32_t *) vnFrameBytesBuf, vnFrameBuffer, nPixelsPerFrame);
            } else {
               letohl((uint32_t *) vnFrameBytesBuf, vnFrameBuffer, nPixelsPerFrame);
            }
            break;
            
         default:
            fprintf(stderr, "*** bin_frames_time/ReadFrame: Invalid pixel size.\n");
            return -1;
      }
      return 0;
   }
}

void WriteFrame(uint64_t *vnFrameBuffer, FILE *hFile, uint8_t *vnFrameBytesBuf, unsigned long int nPixelsPerFrame, long unsigned int nBytesPerPixel, bool bBigEndian) {
   int nIndex;
   
   /* Convert endian-ness */
   switch (nBytesPerPixel) {
      case 1:
         for (nIndex = 0; nIndex < nPixelsPerFrame; nIndex++) {
            vnFrameBytesBuf[nIndex] = (uint8_t) vnFrameBuffer[nIndex];
         }
         
         break;
         
      case 2:
         if (bBigEndian) {
            htobes(vnFrameBuffer, (uint16_t *) vnFrameBytesBuf, nPixelsPerFrame);
         } else {
            htoles(vnFrameBuffer, (uint16_t *) vnFrameBytesBuf, nPixelsPerFrame);
         }
         break;
         
      case 4:
         if (bBigEndian) {
            htobel(vnFrameBuffer, (uint32_t *) vnFrameBytesBuf, nPixelsPerFrame);
         } else {
            htolel(vnFrameBuffer, (uint32_t *) vnFrameBytesBuf, nPixelsPerFrame);
         }
         break;
      default:
         errprintf("*** bin_frames_time/WriteFrame: Invalid pixel size.\n");
         exit(-1);
   }
   
   /* - Write buffer */
   if (fwrite(vnFrameBytesBuf, 1, nBytesPerPixel * nPixelsPerFrame, hFile) != nBytesPerPixel * nPixelsPerFrame) {
      errprintf("*** bin_frames_time/WriteFrame: Could not write to output file.\n   Reason: [%s]\n", strerror(errno));
      exit(-1);
   }
}

/* CheckLongArg: Try to parse a long integer from a text string, and die if we can't */
long int CheckLongArg(const char *strText) {
   char 		*strTail;
   long int	nNum;
   
   /* - Try to extract the number */
   nNum = strtol(strText, &strTail, 0);
   
   if (strText == strTail) {
      errprintf("*** bin_frames_time: Could not parse argument [%s] for a long int.\n", strText);
      exit(-1);
   }
   
   return nNum;
}


int TemporalBinning(const char *strProgramName, const char *strInputFile, const char *strOutputFile,
        long unsigned int nXPixels, long unsigned int  nYPixels, long unsigned int  nChannels, long unsigned int nBytesPerPixel, long unsigned int  nFramesPerBin,
        bool bBigEndian) {
   FILE					*hInputFile, *hOutputFile;
   unsigned long int	nGlobalFrameIndex, nPixelsPerFrame, nFramesThisBin = 0;
   uint8_t				*vnFrameBytesBuf;
	uint64_t				*vnBinBuf, *vnThisFrame;
   int nSampleIndex;
   

   /* -- Try to open the files */
   
   if ((hInputFile = fopen(strInputFile, "r" FORCE_BINARY_FILE)) == NULL) {
      errprintf("*** %s/TemporalBinning: Could not open file [%s] for reading.\n   Reason: [%s]\n",
                strProgramName, strInputFile, strerror(errno));
      return -1;
   }
   
   if ((hOutputFile = fopen(strOutputFile, "w" FORCE_BINARY_FILE)) == NULL) {
      fclose(hInputFile);
      
      errprintf("*** %s/TemporalBinning: Could not open file [%s] for writing.\n   Reason: [%s]\n",
              strProgramName, strOutputFile, strerror(errno));
      
      return -1;
   }
   
   
   /* -- Allocate the Bufs */
   
   nPixelsPerFrame = nXPixels * nYPixels * nChannels;
   
   if ((vnBinBuf = malloc(sizeof(uint64_t) * nPixelsPerFrame)) == NULL) {
      fclose(hInputFile);
      fclose(hOutputFile);
      
      errprintf("*** %s/TemporalBinning: Could not allocate bin buffer.\n", strProgramName);
      return -1;
   }
   
   /* - Clear the buffer */
   memset(vnBinBuf, 0, sizeof(uint64_t) * nPixelsPerFrame);
   
   if ((vnThisFrame = (uint64_t *) malloc(sizeof(uint64_t) * nPixelsPerFrame)) == NULL) {
      fclose(hInputFile);
      fclose(hOutputFile);
      
      errprintf("*** %s/TemporalBinning: Could not allocate frame buffer.\n", strProgramName);
      exit(-1);
   }
   
   if ((vnFrameBytesBuf = (uint8_t *) malloc(sizeof(uint8_t) * nBytesPerPixel * nPixelsPerFrame)) == NULL) {
      fclose(hInputFile);
      fclose(hOutputFile);
      free(vnThisFrame);
      
      errprintf("*** %s/TemporalBinning: Could not allocate frame buffer.\n", strProgramName);
      return -1;
   }
   
   
   /* -- Bin frames, loop over frames */
   nGlobalFrameIndex = 0;
   while (!ReadFrame(hInputFile, vnThisFrame, vnFrameBytesBuf, nPixelsPerFrame, nBytesPerPixel, bBigEndian)) {
      /* - Accumulate frame */
      for (nSampleIndex = 0; nSampleIndex < nPixelsPerFrame; nSampleIndex++) {
         vnBinBuf[nSampleIndex] += (uint64_t) vnThisFrame[nSampleIndex];
      }
      
      nFramesThisBin++;
      
      if (nFramesThisBin == nFramesPerBin) {
         /* - Normalise frame */
         for (nSampleIndex = 0; nSampleIndex < nPixelsPerFrame; nSampleIndex++) {
            vnBinBuf[nSampleIndex] /= nFramesThisBin;
         }
         
         WriteFrame(vnBinBuf, hOutputFile, vnFrameBytesBuf, nPixelsPerFrame, nBytesPerPixel, bBigEndian);
         
         nFramesThisBin = 0;
         memset(vnBinBuf, 0, sizeof(uint64_t) * nPixelsPerFrame);

         nGlobalFrameIndex++;
      }
   }

    printf("--- %s: Wrote [%ld] frames.\n", strProgramName, nGlobalFrameIndex);
   
   
   /* -- Free memory */
   
   free(vnBinBuf);
   free(vnThisFrame);
   free(vnFrameBytesBuf);
   
   
   /* -- Close files */
   
   fclose(hInputFile);
   fclose(hOutputFile);
   
   return 0;
}


/* -- main */

int main(int nArgc, char *strArgv[]) {
   /* -- Local parameters and arguments */
   const char 				*strProgramName, *strInputFile, *strOutputFile;
   unsigned long int		nXPixels, nYPixels, nChannels, nBytesPerPixel, nFramesPerBin;
   uint8_t	bBigEndian;
   
   /* -- Check number of arguments */
   if (nArgc < 8) {
      printf("*** %s: Incorrect usage.\n\n", strArgv[0]);
      printf("%s - Average frames in a linear binary file with temporal binning\n", strArgv[0]);
      printf("Usage: %s 'strInputFile' 'strOutputFile' 'nXPixels' 'nYPixels' 'nChannels' 'nBytesPerPixel' 'nFramesPerBin' 'bBigEndian'\n\n", strArgv[0]);
      printf("   'strInputFile' is the file to process.  It is composed of a number of frames,\n");
      printf("   each composed of ['nXPixels' 'nYPixels' 'nChannels'] pixels.  Each pixel is\n");
      printf("   'nBytesPerPixel' is size.  This file will be binned by averaging corresponding\n");
      printf("   pixels for 'nFramesPerBin' frames each, and written to 'strOutputFile' in the\n");
      printf("   same pixel format as the input file.\n\n");
      
      /* - Return an error code */
      exit(-1);
   }
   
   /* -- Parse arguments */
   strProgramName = strArgv[0];
   strInputFile = strArgv[1];
   strOutputFile = strArgv[2];
   nXPixels = CheckLongArg(strArgv[3]);
   nYPixels = CheckLongArg(strArgv[4]);
   nChannels = CheckLongArg(strArgv[5]);
   nBytesPerPixel = CheckLongArg(strArgv[6]);
   nFramesPerBin = CheckLongArg(strArgv[7]);
   bBigEndian = (bool) CheckLongArg(strArgv[8]);
   
/*
   printf(	"Called with [%s] [%s] [%s] [%ld] [%ld] [%ld] [%ld] [%ld] [%d]\n",
           strProgramName, strInputFile, strOutputFile, nXPixels, nYPixels, nChannels,
           nBytesPerPixel, nFramesPerBin, bBigEndian);
*/
   
   /* -- Perform the binning */
   exit(TemporalBinning(	strProgramName, strInputFile, strOutputFile,
                           nXPixels, nYPixels, nChannels, nBytesPerPixel, nFramesPerBin,
                           bBigEndian));
}


/* -- Matlab MEX interfaction function */

#ifdef MATLAB_MEX_FILE

void mexFunction(int nLHS, mxArray *pLHS[], int nRHS, const mxArray *pRHS[]) {
   /* -- Local variables */
   char						*strInputFile, *strOutputFile;
   unsigned long int		nXPixels, nYPixels, nChannels, nBytesPerPixel, nFramesPerBin;
   bool						bBigEndian;
   int						nLength;
   
   /* -- Check arguments */

   if (nRHS < 8) {
      mexErrMsgIdAndTxt("FocusStack:bin_frames_time:InvalidArgument",
              "*** bin_frames_time: Error: Incorrect usage.");
   }
   
   if (!mxIsChar(pRHS[0])) {
      mexErrMsgIdAndTxt("FocusStack:bin_frames_time:InvalidArgument",
              "*** bin_frames_time: Error: 'strInputFile' must be a string.");
   }
   
   if (!mxIsChar(pRHS[1])) {
      mexErrMsgIdAndTxt("FocusStack:bin_frames_time:InvalidArgument",
              "*** bin_frames_time: Error: 'strOutputFile' must be a string.");
   }
   
   if (!mxIsDouble(pRHS[2]) || !mxIsScalar(pRHS[2]) || mxIsComplex(pRHS[2]) ||
           !mxIsDouble(pRHS[3]) || !mxIsScalar(pRHS[3]) || mxIsComplex(pRHS[3]) ||
           !mxIsDouble(pRHS[4]) || !mxIsScalar(pRHS[4]) || mxIsComplex(pRHS[4]) ||
           !mxIsDouble(pRHS[5]) || !mxIsScalar(pRHS[5]) || mxIsComplex(pRHS[5]) ||
           !mxIsDouble(pRHS[6]) || !mxIsScalar(pRHS[6]) || mxIsComplex(pRHS[6]) ||
           !mxIsScalar(pRHS[7]) || mxIsComplex(pRHS[7])) {
      mexErrMsgIdAndTxt("FocusStack:bin_frames_time:InvalidArgument",
              "*** bin_frames_time: Error: 'nXPixels, 'nYPixels', 'nChannels', 'nBytesPerPixel' and 'nFramesPerBin' must be scalar real doubles; 'bBigEndian' must be a real scalar.");
   }
   
   
   /* -- Parse arguments */
   
   nLength = mxGetN(pRHS[0])*sizeof(mxChar)+1;
   strInputFile = malloc(nLength);
   mxGetString(pRHS[0], strInputFile, nLength);
   
   nLength = mxGetN(pRHS[1])*sizeof(mxChar)+1;
   strOutputFile = malloc(nLength);
   mxGetString(pRHS[1], strOutputFile, nLength);
   
   nXPixels = (unsigned long int) mxGetScalar(pRHS[2]);
   nYPixels = (unsigned long int) mxGetScalar(pRHS[3]);
   nChannels = (unsigned long int) mxGetScalar(pRHS[4]);
   nBytesPerPixel = (unsigned long int) mxGetScalar(pRHS[5]);
   nFramesPerBin = (unsigned long int) mxGetScalar(pRHS[6]);
   bBigEndian = (bool) mxGetScalar(pRHS[7]);
   
/*
   printf(	"Called with [%s] [%s] [%ld] [%ld] [%ld] [%ld] [%ld] [%d]\n",
           strInputFile, strOutputFile, nXPixels, nYPixels, nChannels,
           nBytesPerPixel, nFramesPerBin, bBigEndian);
*/
   

   /* -- Call processor function */
   
   TemporalBinning("bin_frames_time", strInputFile, strOutputFile,
           nXPixels, nYPixels, nChannels, nBytesPerPixel, nFramesPerBin,
           bBigEndian);
   
   
   /* -- Free arguments */
   
   free(strInputFile);
   free(strOutputFile);
}

#endif



/* --- END of bin_frames_time.c --- */
