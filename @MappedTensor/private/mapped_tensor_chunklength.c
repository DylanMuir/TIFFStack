/* mapped_tensor_chunklength.c - MEX FUNCTION Calcluate the length of a chunk
 *
 * Usage: [nChunkLength] = mapped_tensor_chunklength(vfDiffs, nStartingIndex)
 *
 * Author: Dylan Muir <muir@hifo.uzh.ch>
 * Created: 17th July, 2012
 */

/* - Fix char16_t definition bug in OS X 10.9 */
#ifndef char16_t
	#include <stdint.h>
	typedef uint16_t char16_t;
#endif
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   /* - Local variables */
   double   *vfDiffs, fComparison, *pfChunkLength;
   mwSize   uNumelDiffs, uStartingIndex,
            uIndex;
   
   /* - Get arguments and sizes */
   uNumelDiffs = mxGetM(prhs[0]) * mxGetN(prhs[0]);
   vfDiffs = mxGetPr(prhs[0]);
   uStartingIndex = (mwSize) *mxGetPr(prhs[1]);   
   
   /* - Allocate destination value */
   if ((plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL)) == NULL) {
      mexErrMsgIdAndTxt("MappedTensor:mapped_tensor_chunklength:Memory",
                        "Could not allocate memory.");
   }
   pfChunkLength = mxGetPr(plhs[0]);
   
   
   /* - Get comparison value */
   fComparison = vfDiffs[uStartingIndex-1];
   
   /* - Scan through data, find first non-match */
   for (uIndex = uStartingIndex; uIndex < uNumelDiffs; uIndex++) {
      if (fComparison != vfDiffs[uIndex]) {
         *pfChunkLength = (double) (uIndex - uStartingIndex + 1);
         return;
      }
   }
   
   /* - Fell through, so the chunk ends at the end of the array */
   *pfChunkLength = (double) (uNumelDiffs - uStartingIndex + 1);
}

/* --- END of mapped_tensor_repsum.c --- */
