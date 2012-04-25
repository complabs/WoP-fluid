//////////////////////////////////////////////////////////////////////////////////////////
// The Matrix class method evaluator. To be used with regression tester.
//
//  Usage:
//    result =  WoP_Eval( operation, A, B )
//    [result,trace]  =  WoP_Eval( operation, A, B )
//
//  where the 'result' depends of the operation and 'trace' is the character 
//  representation of an execetued C++ statement. Input variables A and B are optional
//  having [] as default values. Operation may be '+', '-', '*', '/', '^' etc.
//
//  Filename: WoP_Eval.cpp
//  Revision: 0.3
//  Date:     2012-04-12
//  Author:   Mikica B Kocic 

#include "mex.h"
#include "WoP_Matrix.h"

#include <cstring> // for strcasecmp/stricmp

#if defined(_MSC_VER)
    #define strcasecmp stricmp
#endif

// Compares strings disregarding case of alphabetic characters
inline bool isequal( const char* str1, const char* str2 )
{
    return 0 == strcasecmp( str1, str2 );
}

//////////////////////////////////////////////////////////////////////////////////////////
// The main entry routine for the MATLAB function WoP_Solver.
//
void mexFunction
(
    int nArgOut, mxArray* argOut[], // left-hand side (varargout) of the function
    int nArgIn, const mxArray* argIn[]  // right-hand side (varargin) of the function
    )
{
    if ( nArgIn <= 1 ) {
        mexWarnMsgTxt( "Usage: result = WoP_TestSanity( operation, A, B )" );
        return;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Parse arguments
    
    char opId[ 256 ] = { 0 }; // The first argument: operation ID; not optional
    {
        const mxArray* arg = argIn[0];
        if ( mxIsChar( arg ) ) {
            mxGetString( arg, opId, sizeof( opId ) - 1 );
        } else {
            mexErrMsgTxt( "Operation must be a string." );
        }
    }

    HeapStack stack;
    // HeapStack::ShowStatistics ();

    Matrix A( nArgIn, argIn, /*arg#*/ 1, "A" ); // 2nd argument, if any
    Matrix B( nArgIn, argIn, /*arg#*/ 2, "B" ); // 3rd argument, if any
    Matrix C( nArgIn, argIn, /*arg#*/ 3, "C" ); // 4th argument, if any
    Matrix D( nArgIn, argIn, /*arg#*/ 4, "D" ); // 5th argument, if any
    Matrix Result;

    //////////////////////////////////////////////////////////////////////////////////////
    // Execute operation
    
    // Keep track of what we are doing...
    const char* title = "";
    #define trace(x) ( title = #x, (x) )

    #define testOpUnary(opName) \
        if ( isequal( opId, #opName ) ) { Result = trace( A.opName () ); }

    #define testOpUnaryScalar(opName) \
        if ( isequal( opId, #opName ) ) { Result = Matrix( 1, 1, trace( A.opName () ) ); }

    #define testOpBinary(opName) \
        if ( isequal( opId, #opName ) ) { Result = trace( A.opName( B ) ); }
    
    if ( isequal( opId, "+" ) ) 
    {
        Result = trace( A + B );
    } 
    else if ( isequal( opId, "-" ) )
    {
        if ( ! B.IsEmpty () ) {
            Result = trace( A - B );
        } else {
            Result = trace( -A );
        }
    }
    else if ( isequal( opId, "*" ) )
    {
        Result = trace( A * B );
    }
    else if ( isequal( opId, "/" ) )
    {
        Result = trace( A / B );
    }
    else if ( isequal( opId, "^" ) )
    {
        Result = trace( A ^ B ); // A ^ B == pow(A,B)
    }
    else testOpBinary( Cross )
    else testOpBinary( CrossCross )
    else testOpBinary( DotProd_ByRow )
    else testOpUnaryScalar( SquaredNorm )
    else testOpUnaryScalar( Norm )
    else testOpUnaryScalar( Sum )
    else testOpUnary( SquaredNorm_ByRow )
    else testOpUnary( Norm_ByRow )
    else testOpUnary( Sum_ByRow )
    else testOpUnary( Sum_ByColumn )
    else if ( isequal( opId, "Row" ) && B.IsScalar () )
    {
        int b = int(B(0)) - 1;
        Result = trace( A.Row( b ) );
    }
    else if ( isequal( opId, "SetRow" ) && B.IsScalar () )
    {
        int b = int(B(0)) - 1;
        Result = trace( A.SetRow( b, C ) );
    }
    else if ( isequal( opId, "Column" ) && B.IsScalar () )
    {
        int b = int(B(0)) - 1;
        Result = trace( A.Column( b ) );
    }
    else if ( isequal( opId, "SetColumn" ) && B.IsScalar () )
    {
        int b = int(B(0)) - 1;
        Result = trace( A.SetColumn( b, C ) );
    }
    else if ( isequal( opId, "Interp" ) )
    {
        // Returns interpolation of heightfield A with x = B(0) and y = C(0)
        // and HF_dim = D
        Result = trace( Matrix( 1, 1, A.Interpolate( B(0), C(0), D ) ) );
    }
    else
    {
        Matrix::SevereError( "Wop:TestSanity", "Unrecognized operation '%s'.", opId );
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Return results

    if ( nArgOut == 0 ) {
        if ( nArgIn >= 2 ) disp( A );
        if ( nArgIn >= 3 ) disp( B );
        if ( nArgIn >= 4 ) disp( C );
        if ( nArgIn >= 5 ) disp( D );
        Result.Display( title );
    }
    if ( nArgOut >= 1 ) {
        argOut[0] = Result;
    }
    if ( nArgOut >= 2 ) {
        argOut[1] = mxCreateString( title );
    }
}