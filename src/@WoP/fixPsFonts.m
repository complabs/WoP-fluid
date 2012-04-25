%% Fixes Postscript fonts
% Fixes MATLAB's choice of base PS fonts, after using print command to export
% figures to PS files; replaces 'Helvetica' with 'Arial' and 'Times-Roman' 
% with 'Times New Roman'.
%
% See http://www.mathworks.se/help/techdoc/creating_plots/f3-103191.html#f3-96850)
% for font support in PS in MATLAB.
%
% Filename: @WoP/fixPsFonts.m
% Revision: 0.1
% Date:     2012-03-19
% Author:   Mikica B Kocic

function fixPsFonts( filename )

    if length( filename ) < 4 || ~strcmpi( filename(end-3:end), '.eps' )
        filename = [ filename, '.eps' ];
    end

    % Read in the EPS file
    fid = fopen( filename );
    ff = fread( fid, '*char' )';   % ff = char(fread(fid))';
    fclose( fid );   

    % Replace MATLAB fonts used in EPS to Windows
    % See allowed fontnames in %windir%\fonts\AdobeFnt.lst
    %
    conv_list = { ...
        '/Helvetica-BoldOblique',  '/Arial-BoldItalicMT'; ...
        '/Helvetica-Bold',         '/Arial-BoldMT'; ...
        '/Helvetica-Oblique',      '/Arial-ItalicMT'; ...
        '/Helvetica',              '/ArialMT'; ...
        '/Times-BoldItalic',       '/TimesNewRomanPS-BoldItalicMT'; ...
        '/Times-Bold',             '/TimesNewRomanPS-BoldMT'; ...
        '/Times-Italic',           '/TimesNewRomanPS-ItalicMT'; ...
        '/Times-Roman',            '/TimesNewRomanPSMT'; ...
    };

    for k = 1 : size( conv_list, 1 )
        ff = strrep( ff, conv_list{k,1}, conv_list{k,2} );
    end

    % Rerite the file with new contents
    fid = fopen( filename, 'w' );
    fprintf( fid, '%s', ff );
    fclose( fid );

end % function