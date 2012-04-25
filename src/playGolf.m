%% LAB3: Play Golf (GUI)
% Allows the user to interactivelly configure parameters and run the game.
%
%  Filename: playGolf.m
%  Revision: 0.2
%  Date:     2012-03-31
%  Author:   Mikica B Kocic
%
%#ok<*DEFNU,*INUSL,*INUSD>

function varargout = playGolf( varargin )

    % ------------------------------------------------------------------------------------
    % Begin initialization code - DO NOT EDIT
    % ------------------------------------------------------------------------------------
    gui_Singleton = 1;

    gui_State = struct( 'gui_Name',       mfilename, ...
                        'gui_Singleton',  gui_Singleton, ...
                        'gui_OpeningFcn', @playGolf_OpeningFcn, ...
                        'gui_OutputFcn',  @playGolf_OutputFcn, ...
                        'gui_LayoutFcn',  [] , ...
                        'gui_Callback',   [] );

    if nargin && ischar( varargin{1} )
        gui_State.gui_Callback = str2func( varargin{1} );
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn( gui_State, varargin{:} );
    else
        gui_mainfcn( gui_State, varargin{:} );
    end
    % ------------------------------------------------------------------------------------
    % End initialization code - DO NOT EDIT
    % ------------------------------------------------------------------------------------

% ========================================================================================
%% Nested functions ----------------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% » initialize_GUI( hObject, handles )
% Initialzes values and strings of the GUI components from simObj properties

function initialize_GUI( hObject, handles )

    handles.simObj.v_fluid  = [ 2 + randn * 6, -2 + randn * 3, 0 ]; 

    set( handles.combo_loft, 'String', handles.params.clubs(:,2) );

    set( handles.text_info, 'String', [  ...
          sprintf( 'Distance:\n\n%10.1f m to hole\n', ...
                   norm( handles.simObj.X - handles.simObj.X_sink ) ), ...
          sprintf( '\nWind:\n  x:%6.1f m/s\n  y:%6.1f m/s', ...
                   handles.simObj.v_fluid(1), handles.simObj.v_fluid(2) ), ...
    ]);
          %sprintf( '  x:%6.1f m\n  y:%6.1f m\n  y:%6.1f m\n', ...
          %         handles.simObj.X(1), handles.simObj.X(2), handles.simObj.X(3) ), ...

    set( handles.slider_direction, 'Value',  handles.params.direction );
    set( handles.slider_speed,     'Value',  handles.params.speed );
    set( handles.combo_loft,       'Value',  handles.params.club );
    set( handles.slider_spin,      'Value',  handles.params.spin );
    set( handles.slider_spin,      'Value',  handles.params.spin );
    set( handles.slider_spin_dir,  'Value',  handles.params.spin_dir );

    set( handles.panel_Stroke, 'Title',  ...
        sprintf( 'Stroke #%d', handles.params.stroke ) ...
    );

    set( handles.panel_direction, 'Title',  ...
        sprintf( 'Direction: %+.0f°', handles.params.direction ) ...
    );

    set( handles.panel_speed, 'Title',  ...
        sprintf( 'Speed: %.0f km/h', handles.params.speed ) ...
    );

    set( handles.panel_spin, 'Title',  ...
        sprintf( 'Spin: %.0f rpm, %+.0f°', ...
            handles.params.spin, handles.params.spin_dir ) ...
    );

    % Setup button Start depending on wheter the WoP simulation is running or not 
    % and if not running on value of t_0
    if handles.simObj.Running
        set( handles.button_Start, 'Enable', 'off', ...
            'String', 'Running...', ...
            'ForegroundColor', [ 0.5, 0.5, 0.5 ], 'FontWeight', 'normal' ...
        );
    elseif handles.params.completed
        set( handles.button_Start, 'Enable', 'on', ...
            'String', 'New Game', ...
            'ForegroundColor', [ 1, 0, 0 ], 'FontWeight', 'bold' ...
        );
    else
        set( handles.button_Start, 'Enable', 'on', ...
            'String', 'Hit!', ...
            'ForegroundColor', [ 0, 0, 0 ], 'FontWeight', 'normal' ...
        );
    end

    % Follow the ball with camera; however, if completed go back to default view

    L = abs([ handles.simObj.BoxWidth, handles.simObj.BoxDepth, ...
              handles.simObj.BoxHeight ]);

    if handles.params.completed

        set( handles.axes_Golf, ...
            'CameraPosition',   [ 0.5 * L(1), -1 * L(2), 1.4* L(3) ], ...
            'CameraTarget',     [ 0.5 * L(1), 0.5 * L(2), 0.3 * L(3) ], ...
            'CameraUpVector',   [ 0, 0, 1 ], ...
            'CameraViewAngle',  20 ...
        );
    
    else
        % Get displacement and distance between the sink and the golf ball ensuring
        % a minimum distance of 10 m.
        displacement = handles.simObj.X_sink - handles.simObj.X(1,:);
        distance = max( eps, norm(displacement) );
        unit_vector = displacement / distance;
        distance = max( 10, distance ); % min distance = 10 m
        displacement = distance * unit_vector;

        % Look behind from the double displacement at height proportional 
        % to box height/depth
        lookFrom = handles.simObj.X(1,:) - displacement;
        lookFrom(3) = max( 4, 1.4 * distance * abs( L(3) / L(2) ) );

        % Look at the half distance between the sink and the ball
        lookAt = ( handles.simObj.X_sink + handles.simObj.X(1,:) ) / 2;
        lookAt(3) = 0.3 * lookFrom(3);

        set( handles.axes_Golf, 'CameraPosition', lookFrom, 'CameraTarget', lookAt, ...
            'CameraViewAngle', 27 - 5 * ( distance / max(L) ) );
    end

    guidata( hObject, handles ); % Update handles structure

% ----------------------------------------------------------------------------------------
%% » parseEditBox( hObject, handles, field, str2value, isValid, errMsg )
% Parses and validates string from editbox and stores value in 'simObj' field (property).
%
% * hObject:   edit box handle
% * handles:   structure with handles and user data (see GUIDATA)
% * field:     name of the field in handles.simObj structure
% * str2Value: function handle that converts string to param value
% * isError:   function handle that validates param value
% * errMsg:    error message displayed if validation fails

function parseEditBox( hObject, handles, field, str2value, isValid, errMsg )

    val = str2value( get( hObject, 'String' ) );

    if isValid( val )
        set( hObject, 'String', val );
        guidata( hObject, handles ) % Update handles structure
        handles.params.(field) = val;
    else % error:
        uicontrol( hObject );
        errordlg( errMsg, 'Play Golf: Validation Error' );
        set( hObject, 'String', handles.simObj.(field) );
    end
    
    guidata( hObject, handles ) % Update handles structure
    
% ========================================================================================
%% App Event Handlers --------------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% » playGolf_OpeningFcn 
% Executes just before playGolf is made visible.
%
% * hObject:    handle to figure
% * eventdata:  reserved - to be defined in a future version of MATLAB
% * handles:    structure with handles and user data (see GUIDATA)
% * varargin:   command line arguments to playGolf (see VARARGIN)

function playGolf_OpeningFcn( hObject, eventdata, handles, varargin )

    % Choose default command line output for playGolf
    handles.output = hObject;

    % Add new invisible axes (in case someone draws to our figure)
    %axes; set( gca, 'Visible', 'off' );

    set( handles.appFigure, 'Toolbar', 'none' );
    cameratoolbar( 'Show' );

    set( handles.appFigure, 'Name', 'Lab3: The Golf Game (v0.2)' );

    % Delayed initialization of GUI components using single-shot timer (note that at
    % this time in OpeningFcn, uicontrol objects are not still created, so changing
    % their properties will result in error).
    % This sets up the initialization when we are invisible so window can get raised 
    % using the 'playGolf' command.

    if strcmpi( get( hObject, 'Visible' ), 'off' )
        start( timer( 'ExecutionMode', 'singleShot', ...
            'TimerFcn', { ...
                @( obj, event, handles ) ...
                    configuration_Callback( ...
                        hObject, handles ...
                    ), ...
                handles ...
            }...
        ) );
    end

% Initialize new simulation

function configuration_Callback( hObject, handles )

    try
        % Delete existing simulation, if any
        clearvars handles.simObj

        handles.params.stroke    = 1; % stroke number#
        handles.params.direction = 0;
        handles.params.speed     = 50 * 3.6;
        handles.params.club      = 8;
        handles.params.spin      = 3600;
        handles.params.spin_dir  = 0;
        handles.params.completed = false;

        % Clubs lookup table:
        % According to: http://en.wikipedia.org/wiki/Iron_%28golf%29 

        handles.params.clubs = ...
        {
            11, '  Iron #0 (11°)'           ;
            14, '  Iron #1 (14°)'           ;
            17, '  Iron #2 (17°)'           ;
            20, '  Iron #3 (20°)'           ;
            23, '  Iron #4 (23°)'           ;
            26, '  Iron #5 (26°)'           ;
            29, '  Iron #6 (29°)'           ;
            33, '  Iron #7 (33°)'           ;
            37, '  Iron #8 (37°)'           ;
            41, '  Iron #9 (41°)'           ;
            45, '  Pitching Wedge (45°)'    ;
            50, '  Gap Wedge (50°)'         ;
            55, '  Sand Wedge (55°)'        ;
            60, '  Lob Wedge (60°)'         ;
            64, '  Ultra Lob Wedge (64°)'   ;
        };

        handles.simObj = WoP( 'golf' );
        handles.simObj.mainAxes = handles.axes_Golf;

        global oobj, oobj = handles.simObj; %#ok<TLEV>

        initialize_GUI( hObject, handles );

        % Update button status to reflect that we are running simulation
        set( handles.button_Start, 'Enable', 'off', ...
            'String', 'Wait...', 'FontWeight', 'bold' ...
        );

        % Start initial simulation (droping ball on the ground, not real-time)
        handles.simObj.InRealTime = false;
        handles.simObj.fps = 1e-3;
        handles.simObj.Silent = true;
        handles.simObj.v_fluid  = [ 0, 0, 0 ]; 
        handles.simObj.Simulate;

        % Restore to real-time mode
        handles.simObj.InRealTime = true; 
        handles.simObj.fps = 25;

        % Update button status to reflect that we are running simulation
        initialize_GUI( hObject, handles );

        % Set minimum figure size (uses obsoleted JavaFrame property, if it exists)
        if usejava( 'swing' )
            % Get current size as the minimum size
            pos = get( handles.appFigure, 'OuterPosition' );
            minimumSize = java.awt.Dimension( pos(3), pos(4) );
            % Set minimim size for the figure's JavaFrame
            try
                warning( 'off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame' );
                jfig = get( handles.appFigure, 'JavaFrame' );
                jfig.fFigureClient.getWindow.setMinimumSize( minimumSize );
            catch %#ok<CTCH>
            end
        end

    catch err
        disp( getReport( err ) )
    end
    
% ----------------------------------------------------------------------------------------
%% » playGolf_OutputFcn
% Outputs from this function are returned to the command line.
%
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function varargout = playGolf_OutputFcn( hObject, eventdata, handles ) 

    % Get default command line output from handles structure
    if isfield( handles, 'output' )
        varargout{1} = handles.output;
    end

% ========================================================================================
%% Button Event Handlers -----------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% » button_Start Callback

function button_Start_Callback( hObject, eventdata, handles )

    if handles.params.completed
        % New game
        configuration_Callback( hObject, handles );
        return;
    elseif handles.simObj.Running
        % Stop any running simulations
        handles.simObj.StopSimulation ();
        return;
    elseif handles.simObj.t_f < handles.simObj.t_0
        % Nothing to simulate
        return;
    end

    fprintf( '==================================================================\n' );
    fprintf( 'STROKE #%d\n', handles.params.stroke );

    % Get club's loft:
    loft = handles.params.clubs{handles.params.club,1};

    fprintf( 'Using %s\n%8.1f°\n', ...
        strtrim( handles.params.clubs{handles.params.club,2} ), loft );

    % Convert units

    loft  = loft * pi/180;                      % 180 deg = pi rad
    v0    = handles.params.speed     / 3.6;     %   1 m/s = 3.6 km/h
    dir   = handles.params.direction * pi/180;  % 180 deg = pi rad
    w0    = handles.params.spin      * 2*pi/60; %  60 rpm = 2*pi rad/s 
    w0ang = handles.params.spin_dir  * pi/180;  % 180 deg = pi rad

    fprintf( 'Direction:\n%8.1f° (ref to ball)\n', dir * 180/pi );

    % Calc absolute direction relative to 'direction 0 = golf ball to sink direction'
    displacement = handles.simObj.X_sink - handles.simObj.X(1,:);
    if displacement(2) >= 0
        dir = dir + atan( displacement(1) / displacement(2) );
    else
        disp( dir );
        disp( atan( displacement(1) / displacement(2) ) );
        dir = dir + pi + atan( displacement(1) / displacement(2) );
    end

    % Round down near-zero values to 0
    w0ang( abs(w0ang) < 1e-3 ) = 0;
    dir( abs(dir) < 1e-3 ) = 0;

    % Setup initial state variables

    handles.simObj.R = handles.simObj.r_p;
    
    handles.simObj.V = [ v0 * cos(loft) * sin(dir) ,  ... % Vx
                         v0 * cos(loft) * cos(dir) ,  ... % Vy
                         v0 * sin(loft)            ];     % Vz

    handles.simObj.W = [  w0 * cos(w0ang) * cos(dir) ,  ... % Wx
                         -w0 * cos(w0ang) * sin(dir) ,  ... % Wy
                         -w0 * sin(w0ang)            ];     % Wz

    % Round down near-zero vector components to 0
    handles.simObj.V( abs(handles.simObj.V) < 1e-4 ) = 0;
    handles.simObj.W( abs(handles.simObj.W) < 1e-4 ) = 0;
    
    % Advance max. final time
    handles.simObj.t_f = handles.simObj.t_f + 20;  

    set( handles.text_info, 'String', ...
        [ sprintf( 'v_0:\n' ), sprintf( '%8.1f m/s\n',   handles.simObj.V ), ...
          sprintf( 'w_0:\n' ), sprintf( '%8.1f rad/s\n', handles.simObj.W ) ] );

    fprintf( '%8.1f° (ref to north)\n', dir * 180/pi );
    disp( get( handles.text_info, 'String' ) );
      
    % Update button status to reflect that we are running application
    set( handles.button_Start, 'Enable', 'off', ...
        'String', 'Running...', ...
        'FontWeight', 'bold', 'ForegroundColor', [ 0.8, 0, 0 ] ...
    );

    guidata( hObject, handles ) % Update handles structure 
    drawnow; % and refresh window before proceeding with simulation

    % Start simulation
    try
        handles.simObj.Silent = true;
        fprintf( '\n>>> Running: %s\n', handles.simObj.Title );
        Simulate( handles.simObj );
        fprintf( 'Elapsed time: %g s\n\n', handles.simObj.t_real );
    catch err
        % Change global flag indicating that the simulation is not running
        handles.simObj.StopSimulation ();

        % Update current gui from the simObj result data
        drawnow; set( 0, 'CurrentFigure', handles.appFigure );
        initialize_GUI( handles.appFigure, handles );

        % Report recognized exceptions separately
        if strcmp( err.identifier, 'MATLAB:hg:dt_conv:Matrix_to_HObject:BadHandle' ) ...
        || strcmp( err.identifier, 'WoP:feim' )
            disp( '>>> Simulation interrupted' );
            disp([ '>>> ', err.message, ':' ]);
            if ~isempty( err.stack )
                disp( err.stack(1) );
            end
            return
        end

        % If not recognized, reissue the exception
        rethrow( err );
    end

    % Reset ball radius from NaN to 
    handles.simObj.R = handles.simObj.r_p;

    % Check for collisions with box surfaces or a entering the sink
    done_info = [];
    if handles.simObj.track.colc(end,2)
        done_info = 'Out of Bounds!';
        handles.params.completed = true;
        text_info = findall( handles.appFigure, 'Tag', 'WoP_textInfo' );
        if ishandle( text_info )
            set( text_info, 'String', done_info, 'Color', 'r', 'FontWeight', 'bold' );
        end
    elseif handles.simObj.track.colc(end,3)
        done_info = 'In the Hole!';
        handles.params.completed = true;
        text_info = findall( handles.appFigure, 'Tag', 'WoP_textInfo' );
        if ishandle( text_info )
            set( text_info, 'String', done_info, 'Color', 'b', 'FontWeight', 'bold' );
        end
    else
        % All ok; increment stroke number
        handles.params.stroke = handles.params.stroke + 1;
    end

    % Update current gui from the simObj result data
    drawnow; 
    figure( handles.appFigure );
    initialize_GUI( handles.appFigure, handles );

    % If completed, issue done_info message
    if ~isempty( done_info )
        msgbox( done_info, 'Lab3: The Golf Game', 'warn' )
    end

% ========================================================================================
%% EditBox Event Handlers  ---------------------------------------------------------------
% Perform parsing and validation of input data and the update of the respective fields
% in 'simObj' structure.

% ----------------------------------------------------------------------------------------
function edit_speed_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'speed', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Speed must be a non-negative real number.' ... error if not
    );

function edit_spin_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'spin', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) ), ... is valid?
        'Spin must be a real number.' ... error if not
    );

% ========================================================================================
%% CheckBox Event Handlers ---------------------------------------------------------------
% Update respective fields in simObj structure with new checkbox values.

% ----------------------------------------------------------------------------------------
function cb_ManualNp_Callback( hObject, eventdata, handles )

    handles.simObj.ManualNp = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
%% ComboBox Event Handlers ---------------------------------------------------------------

function combo_loft_Callback( hObject, eventdata, handles )

    handles.params.club = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
%% CreateFcn Event Handler (for all objects) ---------------------------------------------
% Executes during object creation, after setting all properties.

function edit_CreateFcn( hObject, eventdata, handles )

    if ispc && isequal( get( hObject, 'BackgroundColor' ), ...
                        get( 0, 'defaultUicontrolBackgroundColor' ) )
    % then            
        set( hObject, 'BackgroundColor', 'white' );
    end

function slider_direction_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function slider_speed_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function combo_loft_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function slider_spin_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function slider_spin_dir_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );


% ========================================================================================
% slider_direction_Callback: Executes on slider movement.
% hObject    handle to slider_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function slider_direction_Callback(hObject, eventdata, handles)

    handles.params.direction = get( handles.slider_direction, 'Value' );
    
    set( handles.panel_direction, 'Title',  ...
        sprintf( 'Direction: %+.0f°', handles.params.direction ) ...
    );
    
    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
% slider_speed_Callback: Executes on slider movement.
% hObject    handle to slider_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function slider_speed_Callback(hObject, eventdata, handles)

    handles.params.speed = get( handles.slider_speed, 'Value' );
    
    set( handles.panel_speed, 'Title',  ...
        sprintf( 'Speed: %.0f km/h', handles.params.speed ) ...
    );

    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
% slider_spin_Callback: Executes on slider movement.
% hObject    handle to slider_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function slider_spin_Callback(hObject, eventdata, handles)

    handles.params.spin = get( handles.slider_spin, 'Value' );
    
    set( handles.panel_spin, 'Title',  ...
        sprintf( 'Spin: %.0f rpm, %+.0f°', ...
            handles.params.spin, handles.params.spin_dir ) ...
    );
    
    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
% slider_spin_dir_Callback: Executes on slider movement.
% hObject    handle to slider_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function slider_spin_dir_Callback(hObject, eventdata, handles)

    handles.params.spin_dir = get( handles.slider_spin_dir, 'Value' );
    
    set( handles.panel_spin, 'Title',  ...
        sprintf( 'Spin: %.0f rpm, %+.0f°', ...
            handles.params.spin, handles.params.spin_dir ) ...
    );
    
    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
%% appFigure_ResizeFcn: Executes when appFigure is resized.
% hObject    handle to appFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function appFigure_ResizeFcn( hObject, eventdata, handles )

    posMax = get( handles.appFigure, 'Position' );

    margin = 5;

    % Move panel_Controls (retaining size)
    pos = get( handles.panel_Controls, 'Position' );
    pos(1) = 2 * margin; % left position
    pos(2) = posMax(4) - pos(4) - margin; % bottom position
    set( handles.panel_Controls, 'Position', pos );
    
    % Right position of panel_Controls
    axesLeft = pos(1) + pos(3) + 2 * margin;

    % Move button_Start (retaining size)
    pos = get( handles.button_Start, 'Position' );
    pos(1) = axesLeft; % left position
    pos(2) = posMax(4) - pos(4) - 2 * margin; % bottom position
    set( handles.button_Start, 'Position', pos );
    
    % Move axes_Golf to fill area bellow button_Start
    pos = get( handles.button_Start, 'Position' );
    pos(1) = axesLeft; % left position
    pos(2) = margin; % bottom position
    pos(3) = posMax(3) - axesLeft - 2 * margin; % axes width
    pos(4) = posMax(4) - 4 * margin; % axes height
    set( handles.axes_Golf, 'Position', pos );

% ========================================================================================
