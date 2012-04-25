%% Exports a tagged text annotation to a file
%  Saves tagged annotation of the current figure to EPS file and deletes object.
%
%  Filename: @WoP/ExportAnnotation.m
%  Revision: 0.2
%  Date:     2012-03-19
%  Author:   Mikica B Kocic
%
%  Arguments:
%  * figh: figure holding annotation
%  * tag: tag of the annotation text
%  * dtype: output format type (for print function)
%  * filename: filename to save annotation

function ExportAnnotation( figh, tag, dtype, filename ) 

    % Parameter: margins added around the annotation text box (in pixels)
    sz.margin = [ 5, 5; 20, 5 ]; % (1,:) = left, bottom; (2,:) = right, top

    % Get annotaiton object by tag
    info.origh = findall( figh, 'Tag', tag );
    if isempty( info.origh )
        return;
    end

    % Create invisible figure with white background with invisible axes
    invh = figure( 'Visible', 'off', 'Color', 'w', ...
        'Position', [ 100, 100, 100, 100 ] );

    try
        axis( [ 0, 1, 0, 1 ] );
        set( gca, 'Visible', 'off', ...
            'ActivePositionProperty', 'position', ...
            'Units', 'normalized', 'Position', [ 0, 0, 1, 1 ] );

        % Create text holding annotation with the same text and interpeter as
        % the tagged object
        info.texth = text( 0, 0, get( info.origh, 'String' ), ...
            'FontName', 'Arial', 'FontSize', 10, ...
            'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left', ...
            'Interpreter', get( info.origh, 'Interpreter' ), 'Margin', 1 );

        % Get extent of the text and size of the figure
        sz.extent = get( info.texth, 'Extent' );
        sz.fig = get( invh, 'Position' );

        % Convert extent to hold width/height of the text in pixel units
        sz.extent = sz.fig(3:4) .* sz.extent(3:4);

        % Set figure size = text size + margins
        sz.fig(3:4) = sz.extent + sum( sz.margin, 1 ); % add margines 

        % Convert margin units to data space units
        sz.margin = bsxfun( @rdivide, sz.margin, sz.extent );

        % Reconfigure position of the figure both on screen and paper accordingly
        set( invh, 'Position', sz.fig, ...
            'PaperPositionMode', 'manual', 'PaperUnits', 'points', ... 
            'PaperSize', sz.fig(3:4), 'PaperPosition', [ 0, 0, sz.fig(3:4) ] ...
        );

        % Move the annotation text also keeping requested margins
        set( info.texth, 'Units', 'data', ...
            'Position', [ sz.margin(1,1), 1 - sz.margin(2,2) ] ...
        );

        % Render the figure to an output file
        print( invh, '-r300', dtype, filename );

        % Fix PS fonts for EPS files
        if length( dtype ) >= 5 && strcmp( dtype(1:5), '-deps' )
            WoP.fixPsFonts( filename );
        end

        % Finally, delete original object
        delete( info.origh );

    catch err
        disp( getReport( err ) );
    end

    close( invh );

end % function
