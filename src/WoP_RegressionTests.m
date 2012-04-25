%% Regression tests for the WoP Matrix class
% Should be executed after any modification to the Matrix class.
%
%  Filename: WoP_RegressionTests.m
%  Revision: 0.14
%  Date:     2012-03-25
%  Author:   Mikica B Kocic

%#ok<*NASGU>

function WoP_RegressionTests

    clear WoP_Eval  % forces reload of the mexFunction

    make  % Compile Wop_Eval first

    for i = 1 : 100

        % Generate test vectors

        nR = round( 10 + 10 * rand );
        nC = round(  5 +  5 * rand );

        A = rand( nR, nC );  % matrix
        B = rand( nR, nC );  % matrix
        R = rand( 1, nC );   % row vector
        C = rand( nR, 1 );   % column vector
        S = rand;            % scalar
        D = [];              % temp matrix

        RegressionTest
    end

    clear WoP_Eval  % forces display of the Matrix Heap statistics

    fprintf( 'WoP_Matrix regression tests have been successfully completed.\n\n' )
    return

    %=====================================================================================
    %% Performs single regression test

    function RegressionTest

        % Test binary operations ---------------------------------------------------------

        TestBinaryOp( '+', '+',  '@plus'    )
        TestBinaryOp( '-', '-',  '@minus'   )
        TestBinaryOp( '*', '.*', '@times'   )
        TestBinaryOp( '/', './', '@rdivide' )
        TestBinaryOp( '^', '.^', '@power', eps )

        % Test unary operations ----------------------------------------------------------

        AEval( 'WoP_Eval( ''SquaredNorm'', A )', ...
               'sum( sum( A .^ 2 ) )', nR * nC * 4 * eps )

        AEval( 'WoP_Eval( ''Norm'', A )', ...
               'sqrt( sum( sum( A .^ 2 ) ) )', nR * nC * 4 * eps )

        AEval( 'WoP_Eval( ''Sum'', A )', ...
               'sum( sum( A ) )', nR * nC * 3 * eps )

        AEval( 'WoP_Eval( ''SquaredNorm_ByRow'', A )', ...
               'sum( A .^ 2, 2 )', nR * nC * 4 * eps )

        AEval( 'WoP_Eval( ''Norm_ByRow'', A )', ...
               'sqrt( sum( A .^ 2, 2 ) )', nC * 4 * eps )

        AEval( 'WoP_Eval( ''Sum_ByRow'', A )', ...
               'sum( A, 2 )', nC * 2 * eps )

        AEval( 'WoP_Eval( ''Sum_ByColumn'', A )', ...
               'sum( A, 1 )', nR * 2 * eps )

        % Test accessing row vectors -----------------------------------------------------

        B = floor( 1 + size( A, 1 ) * rand );
        R = A(B,:);
        AEval( 'WoP_Eval( ''Row'', A, B )', 'R' )

        B = floor( 1 + size( A, 1 ) * rand );
        R = rand( 1, nC );   % random row vector
        AEval( 'WoP_Eval( ''Row'', WoP_Eval( ''SetRow'', A, B, R ), B )', 'R' )

        % Test accessing column vectors --------------------------------------------------

        B = floor( 1 + size( A, 2 ) * rand );
        C = A(:,B);
        AEval( 'WoP_Eval( ''Column'', A, B )', 'C' )

        B = floor( 1 + size( A, 2 ) * rand );
        C = rand( nR, 1 );   % random column vector
        AEval( 'WoP_Eval( ''Column'', WoP_Eval( ''SetColumn'', A, B, C ), B )', 'C' )

        % Test cross products ------------------------------------------------------------

        A = rand( 1, 3 );
        B = rand( 1, 3 );

        AEval( 'WoP_Eval( ''Cross'', A, B )', ...
               'cross( A, B )' )

        AEval( 'WoP_Eval( ''Cross'', A, WoP_Eval( ''Cross'', A, B ) )', ...
               'cross( A, cross( A, B ) )' )

        AEval( 'WoP_Eval( ''CrossCross'', A, B )', ...
               'cross( A, cross( A, B ) )', 2 * eps )

        A = rand( 1, 2 );
        B = rand( 1, 2 );

        AEval( '[ WoP_Eval( ''CrossCross'', A, B ), 0 ]', ...
               'cross( [A,0], cross( [A,0], [B,0] ) )', 2 * eps )

        A = rand;
        B = rand;

        AEval( '[ WoP_Eval( ''CrossCross'', A, B ), 0, 0 ]', ...
               'cross( [A,0,0], cross( [A,0,0], [B,0,0] ) )', 2 * eps )

        % Test that, if |A| == 1, always holds: A x ( A x B ) = A <A,B> - B

        A = rand( 1, 3 );
        A = A / WoP_Eval( 'Norm_ByRow', A ); % Makes |A| = 1
        B = rand( 1, 3 );

        AEval( 'WoP_Eval( ''CrossCross'', A, B )', ...
               [ 'WoP_Eval( ''-'', WoP_Eval( ''*'', A, ', ...
                 'WoP_Eval( ''DotProd_ByRow'', A, B ) ', ...
                 '), B )' ], 2 * eps )

        % Test bilinear interpolation (x,y) ----------------------------------------------

        resolution = 0.5;
        [ x, y ] = ndgrid( 0:resolution:20, 0:resolution:10 );
        A = x + y;
        B = x(end) * rand; % = X
        C = y(end) * rand; % = Y
        D = [ x(1) x(end); y(1) y(end) ]; %#ok<SETNU>

        AEval( 'WoP_Eval( ''Interp'', A, B, C, D )', ...
               'B + C', 100 * eps )
           
        % Test linear interpolation (x,*) ------------------------------------------------

        [ x, y ] = ndgrid( 0:resolution:20, 0:resolution:0 );
        A = x + y;
        B = x(end) * rand; % = X
        C = 0;             % = Y
        D = [ x(1) x(end); y(1) y(end) ];

        AEval( 'WoP_Eval( ''Interp'', A, B, C, D )', ...
               'B + C', 100 * eps )

        % Test linear interpolation (*,y)-------------------------------------------------

        [ x, y ] = ndgrid( 0:resolution:0, 0:resolution:20 );
        A = x + y;
        B = 0;             % = X
        C = y(end) * rand; % = Y
        D = [ x(1) x(end); y(1) y(end) ];

        AEval( 'WoP_Eval( ''Interp'', A, B, C, D )', ...
               'B + C', 100 * eps )
    end

    %=====================================================================================
    %% Tests WoP binary operation
    %  Arguments:
    %  * op1: WoP_Eval operation name, e.g. '*'
    %  * op2: MATLAB operation mnemonic, e.g. '.*'
    %  * op3: MATLAB operation function handle to be used with bsxfun(), e.g. @times
    %  * maxerr: allower max absolute error (0 requires exact results)

    function TestBinaryOp( op1, op2, op3, maxerr )
        if nargin < 4
            maxerr = 0;
        end
        AEval( ['WoP_Eval( ''', op1, ''', A, B )'], ['A ', op2, ' B'], maxerr )
        AEval( ['WoP_Eval( ''', op1, ''', A, S )'], ['A ', op2, ' S'], maxerr )
        AEval( ['WoP_Eval( ''', op1, ''', S, A )'], ['S ', op2, ' A'], maxerr )
        AEval( ['WoP_Eval( ''', op1, ''', A, R )'], ['bsxfun(', op3, ', A, R )'], maxerr )
        AEval( ['WoP_Eval( ''', op1, ''', R, A )'], ['bsxfun(', op3, ', R, A )'], maxerr )
        AEval( ['WoP_Eval( ''', op1, ''', A, C )'], ['bsxfun(', op3, ', A, C )'], maxerr )
        AEval( ['WoP_Eval( ''', op1, ''', C, A )'], ['bsxfun(', op3, ', C, A )'], maxerr )
    end

    %=====================================================================================
    %% Asserts equality of two evaluated expressions
    %  Arguments
    %  * e1: 1st expression to be evaluated
    %  * e2: 2nd expression to be evaluated
    %  * maxerr: allower max absolute error (0 requires exact results)

    function AEval( e1, e2, maxerr )

        if nargin < 3
            maxerr = 0;
        end

        % Calculate expressions
        r1 = eval( e1 );
        r2 = eval( e2 );

        % Just return on equal results and continue with tests
        if isequal( r1, r2 )
            return
        end

        % We have different results.
        % Check if difference is bellow the limit...
        abserr = 0;
        if isequal( size(r1), size(r2) )
            abserr = abs(max(max( r1 - r2 )));
            if abserr == 0  % ignore +/-0
                return
            elseif maxerr ~= 0 && abserr <= maxerr
                % Absolute error is bellow the allowed limit; report warning
                % fprintf( 'Warning: %s abs err %e <= %e\n', e1, abserr, maxerr );
                return
            end
        end

        % Report error and terminate
        disp( '--------------------------------------------------------------------' );

        disp( 'A = ' ); disp( A ); 
        disp( 'B = ' ); disp( B ); 
        disp( 'C = ' ); disp( C ); 
        disp( 'R = ' ); disp( R ); 
        disp( 'S = ' ); disp( S ); 
        disp([ e1, ' =' ]); disp( r1 );
        disp([ e2, ' =' ]); disp( r2 );

        if abserr ~= 0
            disp( 'Error =' ); disp( r1 - r2 );
            fprintf( 'Absolute Error = %g (allowed %g)\n\n', abserr, maxerr );
        end

        error( 'AssertEval:failed', 'Assertion failed\n\n  %s ~= %s', e1, e2 );
    end

end