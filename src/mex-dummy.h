#ifndef _MEX_H_INC
#define _MEX_H_INC

// DUMMY MATLAB's mex.h

#include <cstring> // needs size_t

typedef void*  mxArray;
typedef size_t mwSize;
typedef bool   mxLogical;
typedef int    mwIndex;

enum { mxREAL };

extern int mxIsNaN( ... );
extern int mxIsInf( ... );
extern int mxIsDouble( ...  );
extern int mxIsComplex( ... );
extern int mxIsLogicalScalar( ... );

extern double mxGetEps ();
extern double mxGetNaN ();
extern double mxGetInf ();

extern int        mxGetNumberOfDimensions( ... );
extern size_t     mxGetNumberOfElements( ... );
extern mwSize*    mxGetDimensions( ... );
extern double*    mxGetPr( ... );
extern mxLogical* mxGetLogicals( ... );
extern mxArray*   mexGetVariable( ... );

extern void mxSetM( ... );
extern void mxSetN( ... );
extern void mxSetData( ... );
extern void mxSetField( mxArray *pm, mwIndex index, const char *fieldname, mxArray *pvalue );

extern mxArray* mxCreateStructArray( ... );
extern mxArray* mxCreateDoubleScalar( ... );
extern mxArray* mxCreateDoubleMatrix( ... );
extern mxArray* mxCreateLogicalScalar( ... );

extern void* mxMalloc( ... );
extern void mxFree( ... );
extern void mexMakeMemoryPersistent( ... );

extern int mexPrintf ( ... );
extern void mexErrMsgIdAndTxt( ... );

#endif