%% Exports a figure to a file
% Saves figure to a file like 'print' function, but, in case of vector output
% formats, makes all graphic objects opaque during printing to avoid that output 
% vector files are rendered with Z-buffer or OpenGL, in which case, vector files 
% will contain uncompressered rasterized graphics.
%
% Filename: @WoP/ExportFigure.m
% Revision: 0.3
% Date:     2012-03-19
% Author:   Mikica B Kocic

function ExportFigure( figh, dtype, filename )

    switch lower( dtype )
        case '-dmeta',  isVector = true;
        case '-deps',   isVector = true;
        case '-depsc',  isVector = true;
        case '-deps2',  isVector = true;
        case '-depsc2', isVector = true;
        case '-dill',   isVector = true;
        case '-dpdf',   isVector = true;
        case '-dsvg',   isVector = true;
        otherwise,      isVector = false;
    end

    if ~isVector
        % Export rasterized output in 300 dpi
        print( figh, '-r300', dtype, filename );
    else
        % When exporting in vector format, make all graphic objects opaque.
        % Find all transparent objects ...
        objh = findobj( figh, '-property', 'FaceAlpha' );
        objh = findobj( objh, '-not', 'FaceAlpha', 1 );
        % ... and remember their alpha channel
        objv = get( objh, 'FaceAlpha' );

        % Make objects opaque
        set( objh, 'FaceAlpha', 1 );

        % Render the output file with the given format
        print( figh, '-r300', dtype, filename );

        % Restore objects' transparencies
        if length(objh) == 1
            set( objh, 'FaceAlpha', objv );
        else
            set( objh, {'FaceAlpha'}, objv );
        end
    end

    % Fix PS fonts for EPS files
    if length( dtype ) >= 5 && strcmp( dtype(1:5), '-deps' )
        WoP.fixPsFonts( filename );
    end

end % function