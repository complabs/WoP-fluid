%% make - Recompiles MEX methods for the WoP class (if necessary)
%
%  Filename: make.m
%  Revision: 0.2
%  Date:     2012-03-31
%  Author:   Mikica B Kocic

function make ()

    mflags = '-O';                % Compile with optimizations
    % mflags = '-g -DDEBUG';      % Compile with debug support
    % mflags = '-v -O';           % Verbose compile

    recompile( mflags, 'WoP_Eval',   'WoP_Eval.cpp',   'WoP_Matrix.h' )
    recompile( mflags, 'WoP_Solver', 'WoP_Solver.cpp', 'WoP_Matrix.h' )

end

%% Recompiles mexfile that is dependend on varargin{} source files.
% The mexfile should be specified without a file extension .mex*.
% The first varargin should be always specified and it should contain a name of the 
% base source file from which mexfile is built, e.g. for mexfile 'test'
% argument varagin{1} should be 'test.cpp' in case of C++.

function recompile( mflags, mexfile, varargin )

    if nargin < 3
        error( 'At least one source file should be specified.' );
    end

    import java.io.*

    rebuild = false;

    bin = File([ mexfile, '.', mexext ]);

    if ~bin.exists ()
        % Turn on recompile flag
        rebuild = true;
    else
        % Check input file by file if it is more recent than the output
        for i = 1:length(varargin)
            src = File( varargin{i} );
            if ~src.exists ()
                % Report missing source files
                error([ 'Source file ', varargin{i}, ' does not exist.' ]);
            elseif src.lastModified() > bin.lastModified ()
                % Source file is more recent, turn on recompile flag
                rebuild = true;
            end
        end
    end

    if rebuild
        cmd = [ 'mex ', mflags, ' ', varargin{1} ];
        disp( cmd )
        eval( cmd )
    end
end
