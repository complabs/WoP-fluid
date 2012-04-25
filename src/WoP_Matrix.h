//////////////////////////////////////////////////////////////////////////////////////////
// The Matrix class used by the World of Particles (WoP) Solver mexFunction.
// Created for The Physics of Virtual Environments (2012) Computer Lab 3.
//
//  Design guidelines for the implementation were:
//
//   * Avoid frequent use of mxMalloc's -- it's usage costs time!
//   * Use inline class methods and don't afraid to be verbosive with a[i++] = b[j++],
//     (instead *pa++ = *pb++) and let the optimzer do its magic and optimize the code.
//
//  The resulting implementation is summarized in two classes: Matrix and MatrixHeap.
//
//  All methods are inline and act more like macros, so the compiler optimizer 
//  can do it's job properlly. All temporary matrices are allocated on a separate
//  heap (allocated once by mxMalloc); the heap is recycled att the begining of 
//  each simualation cycle or function.
//
//  Be careful that some of the (mostly private) inline methods do not check array limits.
//
//  Tested with:
//
//  MATLAB R2010a, 32- & 64-bit
//      MS Visual C++ 2010 32- & 64-bit
//      MS Visual C++ 2008 32-bit
//      MS Visual C++ 6.0  32-bit
//
//  MATLAB R2011a, 64-bit 
//      GNU C++ 4.4.5-8 on Debian 64bit
//
//  For performance comparison with MATLAB, see WOP_Solver implementation.
//
//  Filename: WoP_Matrix.cpp
//  Revision: 0.16
//  Date:     2012-04-22
//  Author:   Mikica B Kocic 

#ifndef _WOP_MATRIX_H_INCLUDED
#define _WOP_MATRIX_H_INCLUDED

#include "mex.h"

#include <cmath>      // uses pow()
#include <cstdarg>    // uses va_list
#include <cstring>    // uses memcpy
#include <cstdio>     // uses vsprintf
#include <algorithm>  // uses std::max

//////////////////////////////////////////////////////////////////////////////////////////
// dbgf() macro

#ifdef DEBUG
    #define dbgf  mexPrintf
#elif defined(_MSC_VER) && _MSC_VER < 1400
    static inline void dbgf(...) {/*nop*/}
#else
    #define dbgf(a,...) /*nop*/
#endif

//////////////////////////////////////////////////////////////////////////////////////////
// Macro for verbose display of expressions involving matrices; functions as 
// a wrapper for Matrix::Display() method.

#define disp( x ) ( (x).Display( #x ) )

//////////////////////////////////////////////////////////////////////////////////////////
// Support for missing std::max() template in MS VS < 2005, e.g. MSVC 6.0

#if defined(_MSC_VER) && _MSC_VER < 1400
namespace std
{
    template <class T> inline const T& max( const T& x, const T& y ) { 
        return x < y ? y : x;
    }
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////
// The Matrix heap class. Implements a heap where Matrix data is allocated from.

class MatrixHeap
{
    enum { MaxParcelCount = 256 };
    
    struct Parcel
    {
        double* area;  // Allocated area using mxMalloc
        size_t  size;  // total area size in bytes
        size_t  free;  // free area size in bytes
        double* next;  // pointer to next free space in allocated area

        bool Initialize( size_t sizeInBytes )
        {
            // We are keeping parcels at double boundary so we need to adjust the 
            // requested size to be a multiple of sizeof( double )
            size_t ndbls = sizeInBytes / sizeof( double ) + 1;
            
            free = size = ndbls * sizeof( double );
            area = next = (double*) mxMalloc( size );
            return area != NULL;
        }

        void Delete ()
        {
            if ( area != NULL ) {
                mxFree( area );
            }
        }

        void Recycle ()
        {
            free = size;
            next = area;
        }

        double* Malloc( size_t sizeInBytes )
        {
            // We are keeping parcels at double boundary so we need to adjust the 
            // requested size to be a multiple of sizeof( double )
            size_t ndbls = sizeInBytes / sizeof( double ) + 1;
            sizeInBytes = ndbls * sizeof( double );

            // Advance next pointer and reduce free space
            double* ptr = next;
            next += ndbls;
            free -= sizeInBytes;
            return ptr;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////
    
    const char* mName;    // Heap name
    Parcel mParcels[ MaxParcelCount ]; // List of parcels
    size_t mDefaultSize;  // Default parcel area size in bytes
    size_t mTotalSize;    // Total size in bytes of all parcel areas
    int mParcelCount;     // Number of parcels
    int mCurrentParcel;   // Current parcel from which Malloc() works
    bool mPersistent;     // Keep allocated space persistent between mexFunction calls

    size_t mMaxAllocSize; // Maximum allocated size at once
    double mTotAllocCount; // Total number of Malloc() calls (in lifetime)
    double mTotAllocSize;  // Total allocated size (in lifetime)
    size_t mAllocCount;    // Number of Malloc() calls since last Recycle()
    size_t mAllocSize;     // Allocated size since last Recycle()

    /////////////////////////////////////////////////////////////////////////////////////

    void SevereError( const char* reason );

    void Expand( size_t minimum_size = 0 )
    {
        // Return a recycled parcel
        if ( ++mCurrentParcel < mParcelCount ) {
            mParcels[ mCurrentParcel ].Recycle ();
            return;
        }

        // otherwise, no more parcels to recycle, we need to allocate a new one...
        
        if ( ++mParcelCount >= MaxParcelCount ) {
            --mCurrentParcel; 
            --mParcelCount;
            SevereError( "Maximum number of parcels reached." );
            // SevereError calls mexError and mexFunction should quit at this point
            return;
        }

        size_t size = std::max( minimum_size, mDefaultSize );
        if ( ! mParcels[ mCurrentParcel ].Initialize( size ) )
        {
            // It is severe error if mxMalloc failed
            --mCurrentParcel; 
            --mParcelCount;
            SevereError( "Could not allocate a new parcel." );
            // SevereError calls mexError and mexFunction should quit at this point
            return;
        }

        mTotalSize += size; // We have a new parcel
        
        // Make the parcel persisent and report expansion
        //
        if ( mPersistent ) 
        {
            mexMakeMemoryPersistent( mParcels[ mCurrentParcel ].area );

            mexPrintf(
                "%s + %u = total %u bytes %s\n"
                "Current parcel# = %d, beg = %p, end = %p\n", 
                mName, size, mTotalSize, mPersistent ? "(persisent)" : "",
                mCurrentParcel, mParcels[ mCurrentParcel ].area,
                mParcels[ mCurrentParcel ].area + size
            );
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////

public:

    /////////////////////////////////////////////////////////////////////////////////////

    MatrixHeap( const char* name, 
            bool persistent = false, size_t defaultParcelSize = 256000 )
    {
        // Heap identifier
        mName = name;

        // Parcel state variables
        mPersistent    = persistent;
        mDefaultSize   = defaultParcelSize;
        mTotalSize     = 0;
        mCurrentParcel = -1;
        mParcelCount   = 0;

        // Statistics
        mMaxAllocSize  = 0;
        mTotAllocCount = 0;
        mTotAllocSize  = 0;
        mAllocCount    = 0;
        mAllocSize     = 0;

        // We must have at least one parcel
        Expand ();  
    }

    /////////////////////////////////////////////////////////////////////////////////////

    ~MatrixHeap ()
    {
        if ( ! mPersistent ) {
            // The MATLAB memory management facility automatically deallocates
            // all of a MEX-file's memory allocations when the MEX-file completes.
            return;
        }

        ShowStatistics ();

        // Free all of our allocated memory explicitly
        //
        for ( int i = 0; i < mParcelCount; ++i ) {
            mParcels[ i ].Delete ();
        }
    }

    void ShowStatistics ()
    {
        mexPrintf( 
            "\n%s Statistics:\n"
            "    Size: %u bytes in %d parcels, top usage: %u bytes\n"
            "    %.3g bytes allocated during lifetime by %.3g alloc calls\n"
            "    %u bytes allocated since last recyle by %u alloc calls\n\n", 
            mName, mTotalSize, mParcelCount, mMaxAllocSize, 
            mTotAllocSize, mTotAllocCount, mAllocSize, mAllocCount
        );
    }

    /////////////////////////////////////////////////////////////////////////////////////

    double* Malloc( size_t size )
    {
        // Get a new parcel if there is no free space in the current parcel
        if ( size > mParcels[ mCurrentParcel ].free ) {
            Expand( size );
        }

        // Keep statistics: number of Malloc() calls and allocated size
        ++mAllocCount;     mAllocSize    += size;
        ++mTotAllocCount;  mTotAllocSize += size;

        if ( mAllocSize > mMaxAllocSize ) {
            mMaxAllocSize = mAllocSize;
        }

        // Allocate requested space from the current parcel
        return mParcels[ mCurrentParcel ].Malloc( size );
    }

    /////////////////////////////////////////////////////////////////////////////////////

    void Recycle ()
    {
        // Reset statistics
        mAllocCount = 0;
        mAllocSize  = 0;

        // Setup initial parcel
        mCurrentParcel = 0;
        mParcels[ mCurrentParcel ].free = mParcels[ mCurrentParcel ].size;
        mParcels[ mCurrentParcel ].next = mParcels[ mCurrentParcel ].area;
    }
};

//////////////////////////////////////////////////////////////////////////////////////
// Stack wrapper for the matrix heaps: 1) remembers current heap in constructor
// 2) sets TempHeap as AllocateOnHeap and 3) restores original heap in destructor.

class HeapStack
{
    MatrixHeap* mOriginal;             // Remembered heap (to be restored on destruct)
    static MatrixHeap* mCurrentHeap;   // Current matrix heap; NULL for mxMalloc's
    static MatrixHeap mTempHeap;       // Matrix heap used for interim results

    friend class Matrix;

public:

    HeapStack( MatrixHeap& heap )
    {
        dbgf( "HeapStack: Current heap %p, replaces %p\n", &heap, mCurrentHeap );

        mOriginal = mCurrentHeap;
        mCurrentHeap = &heap;
    }

    HeapStack ()
    {
        dbgf( "HeapStack: Current heap %p, replaces %p\n", &mTempHeap, mCurrentHeap );

        mOriginal = mCurrentHeap;
        mCurrentHeap = &mTempHeap;
        mTempHeap.Recycle ();
    }

    ~HeapStack ()
    {
        dbgf( "~HeapStack: Current heap %p, replaces %p\n", mOriginal, mCurrentHeap );

        mCurrentHeap = mOriginal;
    }

    static MatrixHeap* Heap ()
    {
        return mCurrentHeap;
    }

    static void UseMxMalloc ()
    {
        mCurrentHeap = NULL;
    }

    static void ShowStatistics ()
    {
        mTempHeap.ShowStatistics ();
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
// The Simple Matrix class
// Follows MATLAB convention for organization of data: MATLAB gives the matrix as rows 
// first, then columns (if you were to traverse the matrix linearly), so (i,j) maps to 
// i + j * rowCount.
//
class Matrix
{
    //////////////////////////////////////////////////////////////////////////////////////
    // Public datatypes

public:
    
    enum AllocType    // Kind of matrix allocation
    {
        MxAllocated,  // Data allocated using mxMalloc
        ExternalData, // Points to external data (e.g. from mxArray)
        OnMatrixHeap, // Data allocated on MatrixHeap
        Temporary     // Data allocated on MatrixHeap as the temporary result
    };

    //////////////////////////////////////////////////////////////////////////////////////
    // Private properties

private:
    
    mutable AllocType mType; // Matrix type

    int     mNRows;   // Number of rows
    int     mNCols;   // Number of columns
    double* mData;   // Matrix data

    //////////////////////////////////////////////////////////////////////////////////////
    // Static private properties

    static bool SuppressMxFree;          // Suppresses mxFree() calls in destructors
    static double eps;                   // Floating-point relative accuracy

    /////////////////////////////////////////////////////////////////////////////////////
    // Private methods

    //////////////////////////////////////////////////////////////////////////////////////
    // Allocates space for data (without initialization) and returns size (in octets) 
    // of the allocated space.
    //
    size_t Allocate( AllocType type, int rowCount, int columnCount )
    {
        mType   = type;
        mNRows  = rowCount;
        mNCols  = columnCount;
        mData   = NULL;

        size_t size = MemorySize ();

        if ( size <= 0 ) {
            // do nothing
        }
        else if ( type != MxAllocated && HeapStack::Heap() != NULL )
        {
            mData = HeapStack::Heap()->Malloc( size );
            if ( mType == ExternalData ) {
                SevereError( "Matrix:Allocate", 
                    "Internal error.\n"
                    "Matrix.Alocate is not allowed for ExternalData [ %d × %d ]\n",
                    mNRows, mNCols 
                );
            }

            dbgf( "%p : >>> NEW [ %d × %d ], mType %d, mData = %p, HEAP %p\n", 
                this, mNRows, mNCols, mType, mData, HeapStack::Heap() );
        } 
        else // Allocate using mxMalloc()
        {
            SuppressMxFree = false; // Enable usage of mxFree()
            mType = MxAllocated;
            mData = (double*)mxMalloc( size ); 
            
            if ( mData == NULL ) {
                SevereError( "Matrix::Allocate", "Failed to allocate memory." );
            }

            dbgf( "%p : >>> NEW [ %d × %d ], mType %d, mData = %p, MXALLOC\n", 
                this, mNRows, mNCols, mType, mData );
        }

        return size;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Makes matrix empty
    //
    void Clear ()
    {
        mNRows = mNCols = 0;
        mType  = ExternalData;
        mData  = NULL;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // Releases allocated space and makes matrix empty.
    //
    void Delete ()
    {
        if ( mData != NULL && mType == MxAllocated && ! SuppressMxFree ) 
        {
            dbgf( "%p : MxFree %p\n", this, mData );
            mxFree( mData );
        }

        Clear ();
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Initializes matrix from mxArray
    //
    void Initialize( const mxArray* arg, const char* argDesc )
    {
        Clear (); // Make matrix empty

        if ( arg == NULL ) {
            return;
        }

        // We expect a scalar, vector or a matrix of real numbers.
        //
        int ndim = mxGetNumberOfDimensions( arg );
        if ( ! mxIsDouble( arg ) || mxIsComplex( arg ) || ndim > 2 )
        {
            SevereError( "WoP:solver:invarg",
                "'%s' must be a scalar, vector or a matrix of real numbers.", 
                argDesc
            );
        }

        // Get the matrix dimensions and the first element of the real data
        const mwSize* dims = mxGetDimensions( arg );
        mNRows = ndim >= 1 ? dims[0] : 1;
        mNCols = ndim >= 2 ? dims[1] : 1;
        mType  = ExternalData;
        mData  = mxGetPr( arg );

        dbgf( "%p : >>> From mxArray [ %d × %d ]\n", this, mNRows, mNCols );
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // Fast alternatives to matrix(ix) and matrix(r,c) without checking dimensions

    double& operator [] ( int ix )
    {
        return mData[ ix ];
    }

    const double& ElementNC( int r, int c ) const
    {
        return mData[ r + c * mNRows ];
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Public static methods and properties

public:
    
    /////////////////////////////////////////////////////////////////////////////////////

    static void SevereError( const char* errorid, const char* errormsg, ... )
    {
        SuppressMxFree = true;
        HeapStack::UseMxMalloc ();

        va_list args;
        va_start( args, errormsg );
        char buffer[ 1024 ];
        vsprintf( buffer, errormsg, args );
        va_end( args );

        mexErrMsgIdAndTxt( errorid, "%s", buffer );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Public methods

public:
    
    //////////////////////////////////////////////////////////////////////////////////////
    // Creates an empty matrix
    //
    Matrix ()
    {
        dbgf( "%p : Matrix (empty)\n", this );

        Clear ();
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Creates an uninitialized temporary matrix on the heap
    
    Matrix( AllocType type, int rowCount, int columnCount )
    {
        dbgf( "%p : Matrix( r %d, c %d )%s\n", this, rowCount, columnCount,
                type == Temporary ? ", Temporary" : "" );

        Allocate( type, rowCount, columnCount );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Copy constructor from a matrix. Gets reference to original matrix data.
    //
    Matrix( const Matrix& copyFrom )
    {
        dbgf( "%p : Matrix( %p ) GET REFERENCE\n", this, &copyFrom );

        // Keep Temporary matrix as temporary, otherwise signal that we have reference 
        // to external (not allocated by us) data. (TODO verify this properly)
        mType  = copyFrom.mType == Temporary ? Temporary : ExternalData;

        mNRows = copyFrom.mNRows;
        mNCols = copyFrom.mNCols;
        mData  = copyFrom.mData;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Copies data from original matrix.
    //
    Matrix& operator = ( const Matrix& copyFrom )
    {
        dbgf( "%p : operator = ( %p ) COPY VALUE\n", this, &copyFrom );

        size_t size = MemorySize ();

        if ( mType != ExternalData && size >= copyFrom.MemorySize () )
        { 
            // Reuse allocated space for new data
            //
            mNRows = copyFrom.mNRows;
            mNCols = copyFrom.mNCols;
        }
        else if ( copyFrom.mType == Temporary ) 
        {
            // Just take over data from temporary matrix
            //
            Delete ();
            size = 0; // Suppress memcpy bellow

            dbgf( "%p : Re-referencing %p\n", this, &copyFrom );
            mType  = OnMatrixHeap;
            mNRows = copyFrom.mNRows;
            mNCols = copyFrom.mNCols;
            mData  = copyFrom.mData;
            copyFrom.mType = ExternalData;
        }
        else 
        {
            // Current space is too small, so the matrix should be reallocated 
            //
            Delete ();
            size = Allocate( OnMatrixHeap, copyFrom.mNRows, copyFrom.mNCols );
        }

        if ( size > 0 ) { // Copy data if any
            memcpy( mData, copyFrom.mData, size );
        }

        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Creates an uninitialized matrix
    //
    Matrix( int rowCount, int columnCount )
    {
        dbgf( "%p : Matrix( r %d, c %d )\n", this, rowCount, rowCount );

        Allocate( OnMatrixHeap, rowCount, columnCount );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Creates matrix from double[] array
    //
    Matrix( const double source [], int rowCount, int columnCount = 1 )
    {
        dbgf( "%p : Matrix( double[] %p, r %d, c %d )\n", 
                this, source, rowCount, columnCount );

        Allocate( OnMatrixHeap, rowCount, columnCount );

        if ( source != NULL ) {
            memcpy( mData, source, MemorySize () );
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Creates a matrix from a scalar value
    //
    Matrix( int rowCount, int columnCount, double value )
    {
        dbgf( "%p : Matrix( r %d, c %d, double %g )\n", 
                this, rowCount, columnCount, value );

        Allocate( OnMatrixHeap, rowCount, columnCount );

        for ( int i = 0; i < Length (); ++i ) {
            mData[ i ] = value;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Constructs matrix from variable from specified workspace
    //
    Matrix( const mxArray* var, const char* varName )
    {
        dbgf( "%p : Matrix( mxArray %p, '%s' )\n", this, var, varName );

        Initialize( var, varName );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Constructs matrix from a mexFunction argument given by argument number
    //
    Matrix( int nArgIn, const mxArray* argIn[], int argNo, const char* argDesc )
    {
        dbgf( "%p : Matrix( argIn %d, '%s' )\n", this, argNo + 1, argDesc );

        if ( argNo < nArgIn )
        {
            Initialize( argIn[ argNo ], argDesc );
        }
        else // argument doesn't exist
        {
            dbgf( "Missing argument %d: %s\n", argNo + 1, argDesc );
            Clear (); // Make matrix empty
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Constructs matrix from variable from specified workspace
    //
    Matrix( const char* varName, const char* workspace = "caller" )
    {
        dbgf( "%p : Matrix( var '%s', workspace '%s' )\n", this, varName, workspace );

        Initialize( mexGetVariable( workspace, varName ), varName );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Destructs an instance releasing all allocated data
    //
    ~Matrix ()
    {
        Delete ();
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Reallocates matrix on mxMalloc() heap
    //
    void MoveToMxHeap ()
    {
        size_t size = MemorySize ();

        if ( mType != MxAllocated && mData != NULL && size > 0 ) {
            const double* oldData = mData;
            Allocate( MxAllocated, mNRows, mNCols );
            memcpy( mData, oldData, size );
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Reallocates matrix on different Matrix heap
    //
    void MoveTo( MatrixHeap& heap )
    {
        size_t size = MemorySize ();

        if ( mData != NULL && size > 0 ) {
            const double* oldData = mData;
            mData = heap.Malloc( size );
            memcpy( mData, oldData, size );
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Resizes matrix (without reshaping or chaning the contents)

    Matrix& Resize( int rowCount, int colCount )
    {
        dbgf( "%p : Resize( %d, %d )\n", this, rowCount, colCount );

        size_t oldSize = MemorySize ();
        size_t newSize = rowCount * colCount * sizeof( double );

        if ( mType != ExternalData && oldSize >= newSize ) { 
            // Just reduce size
            mNRows = rowCount;
            mNCols = colCount;
        }
        else { 
            // Current space is too small and it the matrix should be reallocated
            double* oldData = mData;
            Delete ();
            newSize = Allocate( OnMatrixHeap, rowCount, colCount );
            
            if ( oldData != NULL && mData != NULL && oldSize > 0 ) { // Copy data if any
                memcpy( mData, oldData, oldSize );
            }
        }

        return *this;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // Converts an instance to a MATLAB double matrix (allocated on MATLAB heap)
    //
    operator mxArray* ()
    {
        // Note: The array creation with mxCreateDoubleMatrix initializes
        // the array memory by filling it with zeros. The following code avoids
        // this unnecessary initialization.

        // Create an uninitialized empty array
        mxArray* result = mxCreateDoubleMatrix( 0, 0, mxREAL );
        
        // Set dimensions (note that either of mNRows/mNCols may be 0)
        mxSetM( result, mNRows );  // Set M = row count
        mxSetN( result, mNCols );  // Set N = column count

        // Now, allocate the memory and fill it with our data (if not empty)
        //
        size_t size = MemorySize (); 
        if ( size > 0 )
        {
            double* output = (double*)mxMalloc( size ); // Allocate space for output
            memcpy( output, mData, size );  // Copy our data to output
            mxSetData( result, output );    // Set mxArray data to point to output
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator () ( unsigned char* indices, int destRowCount = -1 ) const
    {
        if ( destRowCount < 0 ) {
            destRowCount = 0;
            for ( int i = 0; i < mNRows; ++i )
            {
                if ( indices[ i ] != 0 ) {
                    ++destRowCount;
                }
            }
        }

        Matrix result( Temporary, destRowCount, mNCols );
        int destRow = 0;

        for ( int i = 0; i < mNRows; ++i )
        {
            if ( indices[ i ] == 0 ) {
                continue;
            }

            // Copy to destination row ( dest, 0 )
            double* to  = &result( destRow++, 0 );
            // vector from row ( i, 0 )
            const double* from = &ElementNC( i, 0 );
            const double* end  = &ElementNC( i, mNCols );
            while( from < end ) { 
                *to = *from; 
                to += result.mNRows; from += mNRows; // traverse rows
            }
        }

        // assert destRow == destRowCount
        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsEmpty: Returns true if the instance is an empty array
    //
    bool IsEmpty () const
    {
        return mNRows * mNCols == 0;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsScalar: Returns true if the instance is a scalar (1x1 matrix)
    //
    bool IsScalar () const
    {
        return mNRows == 1 && mNCols == 1;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsVector: Returns true if the instance is a vector (1xM or Nx1 matrix) or a scalar
    //
    bool IsVector () const
    {
        return ( mNRows == 1 && mNCols > 0 ) || ( mNCols == 1 && mNRows > 0 );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsRowVector: Returns true if the instance is a row vector (1xM matrix)
    //
    bool IsRowVector () const
    {
        return mNRows == 1 && mNCols > 0;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsColumnVector: Returns true if the instance is a column vector (Nx1 matrix)
    //
    bool IsColumnVector () const
    {
        return mNCols == 1 && mNRows > 0;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsFullMatrix: Returns true if the instance is non-empty matrix that is not a vector
    //
    bool IsFullMatrix () const
    {
        return mNCols > 1 && mNRows > 1;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsSize: Checks the sizes of each dimension of the matrix
    //
    bool IsSize( int M, int N ) const
    {
        return M == mNRows && N == mNCols;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsValidRow: Returns true if the specified index points to a valid row
    //
    bool IsValidRow( int r ) const
    {
        return 0 <= r && r < mNRows;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // IsValidColumn: Returns true if the specified index points to a valid column
    //
    bool IsValidColumn( int c ) const
    {
        return 0 <= c && c < mNCols;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Length: Returns the length of the (linear) array of the instance
    //
    int Length () const
    {
        return mNRows * mNCols;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // GetM: Gets number of rows in array
    //
    int GetM () const
    {
        return mNRows;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // GetM: Gets number of columns in array
    //
    int GetN () const
    {
        return mNCols;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // MemorySize: Gets memory size (in bytes) occupied by array
    //
    size_t MemorySize () const
    {
        return mNRows * mNCols * sizeof( double );
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // VerifyDims: Verifies sizes of each matrix dimension and reports severe error
    // if dimensions do not agree.
    //
    void VerifyDims( int M, int N, const char* error_message ) const
    {
        if ( IsSize( M, N ) || ( Length () == 0 && M * N == 0 ) ) {
            return; // OK
        }

        SevereError( 
            "WoP:Matrix:VerifyDims:invsize",
            "%s; Required dimensions [ %d × %d ] have only [ %d × %d ]", 
            error_message, M, N, mNRows, mNCols
        ); 
    }

    //////////////////////////////////////////////////////////////////////////////////////

    double& operator ()( int ix )
    {
        if ( ix < 0 || ix >= Length () ) {
            SevereError( 
                "WoP:Matrix:elem:invdim",
                "Linear index %d outside dimensions [ %d × %d ]", ix, mNRows, mNCols
            ); 
        }

        return mData[ ix ];
    }

    const double& operator ()( int ix ) const
    {
        if ( ix < 0 || ix >= Length () ) {
            SevereError( 
                "WoP:Matrix:elem:invdim",
                "Linear index %d outside dimensions [ %d × %d ]", ix, mNRows, mNCols
            ); 
        }

        return mData[ ix ];
    }

    double& operator ()( int r, int c )
    {
        int ix = r + c * mNRows;
        if ( ix < 0 || ix >= Length () ) {
            SevereError( 
                "WoP:Matrix:elem:invdim",
                "Index (%d,%d) outside dimensions [ %d × %d ]", r, c, mNRows, mNCols
            ); 
        }

        return mData[ ix ];
    }

    const double& operator ()( int r, int c ) const
    {
        int ix = r + c * mNRows;
        if ( ix < 0 || ix >= Length () ) {
            SevereError( 
                "WoP:Matrix:elem:invdim",
                "Index (%d,%d) outside dimensions [ %d × %d ]", r, c, mNRows, mNCols
            ); 
        }

        return mData[ ix ];
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Row(): Gets a row vector from the matrix
    //
    Matrix Row( int i ) const
    {
        Matrix result( Temporary, 1, mNCols );

        const double* ptr = &(*this)( i, 0 ); // Vector in row (i), with limits check
        for ( int j = 0; j < mNCols; ++j ) {
            result.mData[ j ] = *ptr;  ptr += mNRows; // traverse row
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // SetRow(): Replaces a row vector in the matrix
    //
    Matrix& SetRow( int i, const Matrix& vector )
    {
        if ( ! vector.IsVector () || mNCols != vector.Length () ) {
            SevereError( 
                "WoP:Matrix:SetRow:invdim",
                "Row must in be replaced with a vector of the same size.\n"
                "Got dimensions: [ %d × %d ] @ %d <-- [ %d × %d ]", 
                mNRows, mNCols, i, vector.mNRows, vector.mNCols
            );
        }

        double* ptr = &(*this)( i, 0 ); // Vector in row (i), with limits check
        for ( int j = 0; j < mNCols; ++j ) {
            *ptr = vector.mData[ j ];  ptr += mNRows; // traverse row
        }

        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // SetRow(): Sets a row vector in the matrix to scalar velue
    //
    Matrix& SetRow( int i, double value )
    {
        double* ptr = &(*this)( i, 0 ); // Vector in row (i), with limits check
        for ( int j = 0; j < mNCols; ++j ) {
            *ptr = value;  ptr += mNRows; // traverse row
        }

        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Column(): Gets a column vector from the matrix
    //
    Matrix Column( int j ) const
    {
        if ( j < 0 || j >= mNCols ) {
            SevereError( 
                "WoP:Matrix:Column:invdim",
                "Index (%d,%d) outside dimensions [ %d × %d ]", 0, j, mNRows, mNCols
            ); 
        }

        // Return reference to existing column data 
        // (instead of allocating data in a new matrix)
        //
        Matrix result;
        result.mNRows = mNRows;
        result.mNCols = 1;
        if ( mData != NULL ) {
            result.mType = ExternalData;
            result.mData = &mData[ j * mNRows ];
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // SetColumn(): Replaces a column vector in the matrix
    //
    Matrix& SetColumn( int j, const Matrix& vector )
    {
        if ( ! vector.IsVector () || mNRows != vector.Length () ) {
            SevereError( 
                "WoP:Matrix:SetColumn:invdim",
                "Column must be replaced with a vector of the same size."
                "Got dimensions: [ %d × %d ] @ %d <-- [ %d × %d ]", 
                mNRows, mNCols, j, vector.mNRows, vector.mNCols
            );
        }

        double* ptr = &(*this)( 0, j ); // Vector in column (j), with limits check
        for ( int i = 0; i < mNRows; ++i ) {
            *ptr++ = vector.mData[ i ]; // traverse column
        }

        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // SetColumn(): Sets a column vector in the matrix to scalar velue
    //
    Matrix& SetColumn( int j, double value )
    {
        double* ptr = &(*this)( 0, j ); // Vector in column (j), with limits check
        for ( int i = 0; i < mNCols; ++i ) {
            *ptr++ = value; // traverse column
        }

        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator ^ ( double exponent ) const
    {
        if ( exponent == 1 ) {
            return *this;
        }

        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = pow( mData[ i ], exponent );
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix LeftPower( double base ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = pow( base, mData[ i ] );
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix ReplaceZeros( double value ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < Length (); ++i ) {
            result.mData[ i ] = mData[ i ] == 0 ? value : mData[ i ];
        }

        return result;
    }

    Matrix ReplacePositive( double value ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < Length (); ++i ) {
            result.mData[ i ] = mData[ i ] > 0 ? value : mData[ i ];
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Clamps matrix elements to reside between 0 and the right matrix element values.
    // E.g if clamp ==  5, values will be clamped to be in range [0..5]
    // and if clamp == -5, values will be clamped to be in range [-5..0].

    Matrix Clamp( Matrix& clamp ) const
    {
        if ( mNRows != clamp.mNRows && mNCols != clamp.mNCols ) {
            SevereError( "WoP:Matrix:Clamp:invdim", 
                "Number of rows must agree.\n"
                "Got dimensions: [ %d × %d ] / [ %d × %d ]", 
                mNRows, mNCols, clamp.mNRows, clamp.mNCols
            );
        }

        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < Length (); ++i ) 
        {
            double value = clamp.mData[ i ];
            if ( mxIsNaN( mData[ i ] ) ) {
                result.mData[ i ] = value;
            } else if ( value > 0 ) {
                result.mData[ i ] = mData[ i ] > value ? value
                                  : mData[ i ] < 0 ? 0 : mData[ i ];
            } else {
                result.mData[ i ] = mData[ i ] < value ? value
                                  : mData[ i ] > 0 ? 0 : mData[ i ];
            }
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    unsigned char* FindPositive( int& count )
    {
        unsigned char* map = (unsigned char*)HeapStack::mTempHeap.Malloc( 
                sizeof( unsigned char ) * Length () );

        count = 0;
        for ( int i = 0; i < Length (); ++i ) {
            int res = mData[ i ] > 0;
            count += res;
            map[ i ] = res;
        }

        return map;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix& operator = ( double scalar )
    {
        for ( int i = 0; i < Length (); ++i ) {
            mData[ i ] = scalar;
        }

        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator - () const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = - mData[ i ];
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator + ( double scalar ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = mData[ i ] + scalar;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator - ( double scalar ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = mData[ i ] - scalar;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator * ( double scalar ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = mData[ i ] * scalar;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator / ( double scalar ) const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < result.Length (); ++i ) {
            result.mData[ i ] = mData[ i ] / scalar;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    friend inline Matrix operator + ( double scalar, const Matrix& right  );
    friend inline Matrix operator - ( double scalar, const Matrix& right  );
    friend inline Matrix operator * ( double scalar, const Matrix& right  );
    friend inline Matrix operator / ( double scalar, const Matrix& right  );
    
    //////////////////////////////////////////////////////////////////////////////////////
    // BEGIN USAGE of EMIT_BINARY_OP macros
    // How to: define BYNARY_OP(x,y) as e.g. x + y or pow(x,y)

    #define EMIT_BINARY_OP \
    { \
        bool equalRowCount    = mNRows == right.mNRows; \
        bool equalColumnCount = mNCols == right.mNCols; \
        \
        if ( equalRowCount && equalColumnCount ) /* Both have the same dimension */ \
        { \
            /* Just simple linear traversal of all data */ \
            Matrix result( Temporary, mNRows, mNCols ); \
            for ( int i = 0; i < result.Length (); ++i ) { \
                result.mData[ i ] = BINARY_OP( mData[ i ], right.mData[ i ] ); \
            } \
            return result; \
        } \
        else if ( equalRowCount ) \
        { \
            if ( right.IsColumnVector () ) /* Right matrix is a column vector */ \
            { \
                Matrix result( Temporary, mNRows, std::max( mNCols, right.mNCols ) ); \
                double* end = result.mData + result.Length (); \
                double* Z = result.mData; /* Linear traversal of the result */ \
                double* A = mData; /* Linear traversal of the left matrix */ \
                for ( double* B = right.mData; Z < end; B = right.mData ) { \
                    /* B points to mData[0] at the beginning of every row traversal */ \
                    /* and increments with with row traversal */ \
                    for ( int j = 0; j < result.mNRows; ++j ) { \
                        *Z++ = BINARY_OP( *A++, *B++ ); \
                    } \
                } \
                return result; \
            } \
            else if ( IsColumnVector () ) /* Left matrix is a column vector */ \
            { \
                Matrix result( Temporary, mNRows, std::max( mNCols, right.mNCols ) ); \
                double* end = result.mData + result.Length (); \
                double* Z = result.mData; /* Linear traversal of the result */ \
                double* B = right.mData; /* Linear traversal of the right matrix */ \
                for ( double* A = mData; Z < end; A = mData ) { \
                    /* A points to mData[0] at the beginning of every row traversal */ \
                    /* and increments with with row traversal */ \
                    for ( int j = 0; j < result.mNRows; ++j ) { \
                        *Z++ = BINARY_OP( *A++, *B++ ); \
                    } \
                } \
                return result; \
            } \
        } \
        else if ( equalColumnCount ) \
        { \
            if ( right.IsRowVector () ) /* Right matrix is a row vector */ \
            { \
                Matrix result( Temporary, std::max( mNRows, right.mNRows ), mNCols ); \
                double* end = result.mData + result.Length (); \
                double* Z = result.mData; /* Linear traversal of the result */ \
                double* A = mData; /* Linear traversal of the left matrix */ \
                for ( double* B = right.mData; Z < end; ) { \
                    /* B points at current row here */ \
                    for ( int j = 0; j < result.mNRows; ++j ) { \
                        *Z++ = BINARY_OP( *A++, *B ); \
                    } \
                    ++B; /* Increment current B row at the end of every row traversal */ \
                } \
                return result; \
            } \
            else if ( IsRowVector () ) /* Left matrix is a row vector */ \
            { \
                Matrix result( Temporary, std::max( mNRows, right.mNRows ), mNCols ); \
                double* end = result.mData + result.Length (); \
                double* Z = result.mData; /* Linear traversal of the result */ \
                double* B = right.mData; /* Linear traversal of the right matrix */ \
                for ( double* A = mData; Z < end; ) { \
                    /* A points at current row here */ \
                    for ( int j = 0; j < result.mNRows; ++j ) { \
                        *Z++ = BINARY_OP( *A, *B++ ); \
                    } \
                    ++A; /* Increment current A row at the end of every row traversal */ \
                } \
                return result; \
            } \
        } \
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator + ( const Matrix& right ) const
    {
        // Handle case when one of matrices is scalar
        if ( right.IsScalar () ) {
            return *this + right.mData[0];
        } else if ( IsScalar () ) {
            return mData[0] + right;
        }

        #define BINARY_OP(x,y) ((x)+(y))
        EMIT_BINARY_OP
        #undef BINARY_OP

        SevereError( "WoP:Matrix:plus:invdim", 
            "Matrix dimensions must agree.\n"
            "Got dimensions: [ %d × %d ] + [ %d × %d ]", 
            mNRows, mNCols, right.mNRows, right.mNCols
        );
        
        return Matrix ();
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator - ( const Matrix& right ) const
    {
        // Handle case when one of matrices is scalar
        if ( right.IsScalar () ) {
            return *this - right.mData[0];
        } else if ( IsScalar () ) {
            return mData[0] - right;
        }

        #define BINARY_OP(x,y) ((x)-(y))
        EMIT_BINARY_OP
        #undef BINARY_OP
    
        SevereError( "WoP:Matrix:minus:invdim", 
            "Matrix dimensions must agree.\n"
            "Got dimensions: [ %d × %d ] - [ %d × %d ]", 
            mNRows, mNCols, right.mNRows, right.mNCols
        );
        
        return Matrix ();
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator * ( const Matrix& right ) const
    {
        // Handle case when one of matrices is scalar
        if ( right.IsScalar () ) {
            return *this * right.mData[0];
        } else if ( IsScalar () ) {
            return mData[0] * right;
        }

        #define BINARY_OP(x,y) ((x)*(y))
        EMIT_BINARY_OP
        #undef BINARY_OP
    
        SevereError( "WoP:Matrix:multiply:invdim", 
            "Matrix dimensions must agree.\n"
            "Got dimensions: [ %d × %d ] * [ %d × %d ]", 
            mNRows, mNCols, right.mNRows, right.mNCols
        );
        
        return Matrix ();
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator / ( const Matrix& right ) const
    {
        // Handle case when one of matrices is scalar
        if ( right.IsScalar () ) {
            return *this / right.mData[0];
        } else if ( IsScalar () ) {
            return mData[0] / right;
        }

        #define BINARY_OP(x,y) ((x)/(y))
        EMIT_BINARY_OP
        #undef BINARY_OP
    
        SevereError( "WoP:Matrix:rdivide:invdim", 
            "Matrix dimensions must agree.\n"
            "Got dimensions: [ %d × %d ] / [ %d × %d ]", 
            mNRows, mNCols, right.mNRows, right.mNCols
        );
        
        return Matrix ();
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix operator ^ ( const Matrix& right ) const
    {
        // Handle case when one of matrices is scalar
        if ( right.IsScalar () ) {
            return *this ^ right.mData[0];
        } else if ( IsScalar () ) {
            return right.LeftPower( mData[0] );
        }

        #define BINARY_OP(x,y) pow((x),(y))
        EMIT_BINARY_OP
        #undef BINARY_OP

        SevereError( "WoP:Matrix:power:invdim", 
            "Matrix dimensions must agree.\n"
            "Got dimensions: [ %d × %d ] ^ [ %d × %d ]", 
            mNRows, mNCols, right.mNRows, right.mNCols
        );

        return Matrix ();
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // END USAGE of Binary operation macro EMIT_BINARY_OP
    //
    #undef EMIT_BINARY_OP
    
    //////////////////////////////////////////////////////////////////////////////////////
            
    Matrix& operator += ( const Matrix& right )
    {
        *this = *this + right;
        return *this;
    }

    Matrix& operator -= ( const Matrix& right )
    {
        *this = *this - right;
        return *this;
    }

    Matrix& operator *= ( const Matrix& right )
    {
        *this = *this * right;
        return *this;
    }

    Matrix& operator /= ( const Matrix& right )
    {
        *this = *this * right;
        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    double SquaredNorm () const
    {
        double norm = 0;

        for ( int i = 0; i < Length (); ++i ) {
            norm += mData[ i ] * mData[ i ];
        }

        return norm;
    }

    double Norm () const
    {
        return sqrt( SquaredNorm () );
    }

    double Sum () const
    {
        double result = 0;

        for ( int i = 0; i < Length (); ++i ) {
            result += mData[ i ];
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////

    Matrix SquaredNorm_ByRow () const
    {
        Matrix result( Temporary, mNRows, 1 );

        for ( int i = 0; i < mNRows; ++i ) {
            result.mData[ i ] = Row( i ).SquaredNorm ();
        }

        return result;
    }

    Matrix Norm_ByRow () const
    {
        Matrix result( Temporary, mNRows, 1 );

        for ( int i = 0; i < mNRows; ++i ) {
            result.mData[ i ] = sqrt( Row( i ).SquaredNorm () );
        }

        return result;
    }

    Matrix Sum_ByRow () const
    {
        Matrix result( Temporary, mNRows, 1 );

        for ( int i = 0; i < mNRows; ++i )
        {
            double sum = 0;
            const double* p = &ElementNC( i, 0 ); // first element in column
            for ( int j = 0; j < mNCols; ++j ) {
                sum += *p, p += mNRows; // traverse columns
            }
            result.mData[ i ] = sum;
        }

        return result;
    }

    Matrix Sum_ByColumn () const
    {
        Matrix result( Temporary, 1, mNCols );

        for ( int j = 0; j < mNCols; ++j )
        {
            double sum = 0;
            const double* p = &ElementNC( 0, j ); // first element in row
            for ( int i = 0; i < mNRows; ++i ) {
                sum += *p++; // traverse rows
            }
            result.mData[ j ] = sum;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Normalizes row vectors in matrix to be unit vectors.
    // In case when row vectors have zero-norm, uses eps.

    Matrix UnitVector_ByRow () const
    {
        Matrix result( Temporary, mNRows, mNCols );

        for ( int i = 0; i < mNRows; ++i )
        {
            // Calculate the squared norm of the row-vector
            //
            double sum = 0;
            const double* p = &ElementNC( i, 0 ); // first element in column
            for ( int j = 0; j < mNCols; ++j ) {
                sum += (*p) * (*p);
                p += mNRows; // traverse columns
            }

            // Calculate the norm, ensuring that it's not zero
            //
            double norm = sum == 0 ? eps : sqrt( sum );

            // Normalize components i.e. divide components by norm
            //
            p = &ElementNC( i, 0 ); // first element in column (source)
            double* pres = &result( i, 0 ); // first element in column (result)
            for ( int j2 = 0; j2 < mNCols; ++j2 ) {
                *pres = (*p) / norm;
                pres += mNRows, p += mNRows; // traverse columns
            }
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Scalar product: Calculates scalar product of two matrices with row vectors 
    // row by row. Returns column vector.

    Matrix DotProd_ByRow ( const Matrix& right ) const
    {
        if ( mNRows != right.mNRows ) {
            SevereError( "WoP:Matrix:scalarr:invdim", 
                "Number of rows must agree.\n"
                "Got dimensions: [ %d × %d ] / [ %d × %d ]", 
                mNRows, mNCols, right.mNRows, right.mNCols
            );
        }

        Matrix result( Temporary, mNRows, 1 );

        for ( int i = 0; i < mNRows; ++i )
        {
            double sum = 0;
            const double* p1 = &ElementNC( i, 0 ); // first element in column
            const double* p2 = &right.ElementNC( i, 0 );
            for ( int j = 0; j < mNCols; ++j ) {
                sum += *p1 * *p2;
                p1 += mNRows, p2 += mNRows; // traverse columns
            }
            result.mData[ i ] = sum;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Cross(): Calculates cross-product A x B, where A and B are arrays of length 3.

    Matrix Cross( const Matrix& right ) const
    {
        if ( ! IsVector () || ! right.IsVector () 
            || Length () != right.Length () || Length () != 3 )
        {
            SevereError( 
                "WoP:Matrix:Cross:invdim",
                "A and B must vectors with a dimension of length 3."
            );
        }

        Matrix result( Temporary, mNRows, mNCols );
        double *A = mData, *B = right.mData, *C = result.mData;

        C[0] =   A[1] * B[2] - A[2] * B[1];
        C[1] = - A[0] * B[2] + A[2] * B[0];
        C[2] =   A[0] * B[1] - A[1] * B[0];

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // CrossCross(): Calculates double cross-product A x ( A x B ), where
    // A and B are arrays of length between 1 and 3.
    // If the instance is unit vector 'n', then, from the triple vector product: 
    // A x ( B x C ) = B ( A . C ) - C ( A . B )  we get: n x ( n x B ) = 
    // n ( n . B ) - B ( n . n ) = n ( n . B ) - B
    
    Matrix CrossCross( const Matrix& right ) const
    {
        if ( ! IsVector () || ! right.IsVector () 
            || Length () != right.Length () || Length () < 1 || Length () > 3 )
        {
            SevereError( 
                "WoP:Matrix:CrossCross:invdim",
                "A and B must vectors with a dimension of length max 3."
            );
        }

        Matrix result( Temporary, mNRows, mNCols );
        double *A = mData, *B = right.mData, *C = result.mData;

        switch( Length () )
        {
        case 1:
            C[0] = 0;
            break;
        case 2:
            C[0] = - A[1]*A[1]*B[0] + A[0]*A[1]*B[1];
            C[1] =   A[0]*A[1]*B[0] - A[0]*A[0]*B[1];
            break;
        case 3:
            C[0] = - A[1]*A[1]*B[0] - A[2]*A[2]*B[0] + A[0]*A[1]*B[1] + A[0]*A[2]*B[2];
            C[1] =   A[0]*A[1]*B[0] - A[0]*A[0]*B[1] - A[2]*A[2]*B[1] + A[1]*A[2]*B[2];
            C[2] =   A[0]*A[2]*B[0] + A[1]*A[2]*B[1] - A[0]*A[0]*B[2] - A[1]*A[1]*B[2];
            break;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Bilinear interpolation lookup function. 
    // Returns -Inf (instead of NaN) if x, y are out of the height field limits.
    //
    // WARNING: We expect data created by *ndgrid* and not meshgrid!
    //          This means that x relate to rows/ and y to columns.
    //
    double Interpolate( double x, double y, const Matrix& HF_ends
        )
    {
        double error = -mxGetInf (); // Negative Infinity

        if ( IsEmpty () || HF_ends.IsEmpty () ) {
            return error;
        }

        HF_ends.VerifyDims( 2, 2, "Heightfield should be [ xMin xMax ; yMin yMax ]" );
        
        double yMin = HF_ends(1,0), yMax = HF_ends(1,1);
        double xMin = HF_ends(0,0), xMax = HF_ends(0,1);
        
        // HF consits of ( M - 1 ) x ( N - 1 ) equally spaced interpolation segments.
        // Find the lengths of these segments. (set 0 if length could not be determined)
        double xTick = mNRows < 2 ? 0 : ( xMax - xMin ) / ( mNRows - 1 );
        double yTick = mNCols < 2 ? 0 : ( yMax - yMin ) / ( mNCols - 1 );

        if ( xTick == 0 && yTick == 0 ) {
            SevereError( "WoP:Matrix:Interpolate", 
                    "Invalid HF_ends [ %g %g; %g %g ] with HF [ %d × %d ]", 
                    xMin, xMax, yMin, yMax, mNRows, mNCols );
        }

        // Transform y/x into row/column indices
        int i = xTick == 0 ? -1 : int( ( x - xMin ) / xTick );
        int j = yTick == 0 ? -1 : int( ( y - yMin ) / yTick );

        dbgf( 
            "X =%8.4g (%3g..%3g / %3g) -> %3d, "
            "Y =%8.4g (%3g..%3g / %3g) -> %3d\n",
            x, xMin, xMax, xTick, i, y, yMin, yMax, yTick, j
        );

        // Setup cell boundaries
        double x1 = i * xTick, x2 = ( i + 1 ) * xTick;
        double y1 = j * yTick, y2 = ( j + 1 ) * yTick;

        // If hf(x): linear interpolation ------------------------------------------------
        //
        if ( yTick == 0 ) 
        {
            if ( i == mNRows - 1 ) { // at the upper limit
                return ElementNC( i, 0 );
            } else if ( i < 0 || i >= mNRows ) {
                return error;
            }

            return ElementNC( i, 0 ) 
                 + ( ElementNC( i+1, 0 ) - ElementNC( i, 0 ) ) * ( x - x1 ) / xTick;
        }

        // If hf(y): linear interpolation ------------------------------------------------
        //
        if ( xTick == 0 ) 
        {
            if ( j == mNCols - 1 ) { // at the upper limit
                return ElementNC( 0, j );
            } else if ( j < 0 || j >= mNCols ) {
                return error;
            }

            return ElementNC( 0, j ) 
                 + ( ElementNC( 0, j+1 ) - ElementNC( 0, j ) ) * ( y - y1 ) / yTick;
        }

        // If hf(x,y): bilinear interpolation --------------------------------------------
        //
        if ( i == mNRows - 1 && ( j >= 0 && j < mNCols ) ) {
            return ElementNC( i, j ); // at the upper row (y values) limit
        } else if ( j == mNCols - 1 && ( i >= 0 && i < mNRows ) ) {
            return ElementNC( i, j ); // at the upper column (x values) limit
        } else if ( i < 0 || i >= mNRows || j < 0 || j >= mNCols ) {
            return error;
        }

        // Bilinear interpolation:
        // Do linear interpolation first in x, and than in y direction
        // For algorithm see: http://en.wikipedia.org/wiki/Bilinear_interpolation
        //
        double z = ElementNC( i  , j   ) * ( x2 - x ) * ( y2 - y )
                 + ElementNC( i  , j+1 ) * ( x2 - x ) * ( y - y1 )
                 + ElementNC( i+1, j   ) * ( x - x1 ) * ( y2 - y )
                 + ElementNC( i+1, j+1 ) * ( x - x1 ) * ( y - y1 );

        return z / xTick / yTick;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    
    const Matrix& Display( const char* title = NULL ) const
    {
        if ( title != NULL ) {
            mexPrintf( "%s =\n", title );
        }

        if ( IsEmpty () ) {
            mexPrintf( " %7s\n", "[]" );
        }
        else {
            for ( int i = 0; i < mNRows; ++i ) {
                for ( int j = 0; j < mNCols; j++ ) {
                    mexPrintf( " % 7.4e", ElementNC( i, j ) );
                }
                mexPrintf( "\n" );
            }
        }

        if ( title != NULL ) {
            mexPrintf( "\n" );
        }

        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
// Inline operations with a scalar on the left side and a Matrix on the right side

inline Matrix operator + ( double scalar, const Matrix& right  )
{
    Matrix result( Matrix::Temporary, right.mNRows, right.mNCols );

    for ( int i = 0; i < result.Length (); ++i ) {
        result.mData[ i ] = scalar + right.mData[ i ];
    }

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////

inline Matrix operator - ( double scalar, const Matrix& right  )
{
    Matrix result( Matrix::Temporary, right.mNRows, right.mNCols );

    for ( int i = 0; i < result.Length (); ++i ) {
        result.mData[ i ] = scalar - right.mData[ i ];
    }

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////

inline Matrix operator * ( double scalar, const Matrix& right  )
{
    Matrix result( Matrix::Temporary, right.mNRows, right.mNCols );

    for ( int i = 0; i < result.Length (); ++i ) {
        result.mData[ i ] = scalar * right.mData[ i ];
    }

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////

inline Matrix operator / ( double scalar, const Matrix& right  )
{
    Matrix result( Matrix::Temporary, right.mNRows, right.mNCols );

    for ( int i = 0; i < result.Length (); ++i ) {
        result.mData[ i ] = scalar / right.mData[ i ];
    }

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////

inline void MatrixHeap::SevereError( const char* reason )
{
    Matrix::SevereError( "MatrixHeap", reason );
}

//////////////////////////////////////////////////////////////////////////////////////////

inline MatrixHeap& operator <<( MatrixHeap& heap, Matrix& obj )
{
    obj.MoveTo( heap );
    return heap;
}

//////////////////////////////////////////////////////////////////////////////////////////

class Logical
{
    bool value;

public:

    Logical( bool v = false )
    {
        value = v;
    }

    Logical( int nArgIn, const mxArray* argIn[], 
             int argNumber, const char* argDesc, bool defaultValue = false )
    {
        if ( argNumber >= nArgIn )
        {
            // mexPrintf( "Missing argument %d: %s\n", argNo + 1, argId ); 
            value = defaultValue;
            return;
        }

        const mxArray* arg = argIn[ argNumber ];

        if ( arg == NULL ) {
            value = defaultValue;
            return;
        }

        // Parse value from different argument types
        int ndim = mxGetNumberOfDimensions( arg );
        size_t len = mxGetNumberOfElements( arg );

        if ( mxIsLogicalScalar( arg ) && len >= 1 ) 
        {
            value = mxGetLogicals( arg )[0];
            return;
        }
        else if ( mxIsDouble( arg ) && ! mxIsComplex( arg ) || len >= 1 )
        {
            double* data = mxGetPr( arg );
            value = data != NULL && data[0] != 0;
            return;
        }

        // If we couldn't parse logical value, return an error
        Matrix::SevereError( "WoP:solver:parseLogical:invarg",
            "'%s' must be a logical value or a real number.", 
            argDesc
        );
    }
    
    bool operator ! () const
    {
        return !value;
    }

    operator bool () const
    {
        return value;
    }

    mxArray* toMxArray()
    {
        return mxCreateLogicalScalar( value );
    }

    Logical& operator = ( bool v )
    {
        value = v;
        return* this;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
// Static initializers for the HeapStack class

MatrixHeap HeapStack::mTempHeap( "Temp Heap", /*persistent*/ true );

MatrixHeap* HeapStack::mCurrentHeap = NULL;

//////////////////////////////////////////////////////////////////////////////////////////
// Static initializers for the Matrix class

bool Matrix::SuppressMxFree = false;
double Matrix::eps = mxGetEps ();

//////////////////////////////////////////////////////////////////////////////////////////

#endif // _WOP_MATRIX_H_INCLUDED
