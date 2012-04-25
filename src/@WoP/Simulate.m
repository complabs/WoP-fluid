%% The Word of Particles Simulation Engine
% Simulation of a number of spherical objects in a box. Supports both the impulse method
% (instantenous transfer of momenta) or the penalty method (linear springs) as the contact
% model. Implements viscous and Magnus forces. Implements friction with surface contact.
% Implements ground surface as variable height field.
%
%  Usage:
%    Simulate( obj )  or  obj.Simulate ()
%        runs simulation
%
%    track = Simulate( obj )  or   track = obj.Simulate ()
%       runs simulation and returns a structure with the system trajectory in
%       phase-space and (optionally) the trace of the total energy over time
%
%  Glossary:
%   * Jolt:        An instantenous change of the state variable during a time-step.
%   * Projection:  The component of one vector that is in the direction of the other.
%   * Trajectory:  A path in a phase plot that shows how the state of a system
%                 changes over time (a trace of the state variables of a system).
%
%  Filename: @WoP/Simulate.m
%  Revision: 0.3.19
%  Date:     2012-03-28
%  Author:   Mikica B Kocic 
%

function varargout = Simulate( this )

    %%% Timestamp the simulation ---------------------------------------------------------

    this.ID = datestr( now, 'YYYYmmdd_HHMMSS' );

    %%% Announce to the world that the simulation is alive -------------------------------

    this.Running = true; % The simulation will run while this flag is true

    %%% Dump running parameters ----------------------------------------------------------

    if ~this.Silent
        fprintf( '\n---------------------- WoP.Simulation\n\n' );
        disp( this ), disp( 'Enumerations:' ), disp( this.enum );
    end

    %-------------------------------------------------------------------------------------
    %%% Declaration of common variables used across the nested functions
    %
    % All variables used by the WoP_Solver must be declared here,
    % othewise mexFunction will not see them!
    %
    % Background: Using nested functions instead of class methods and declaring common 
    % variables used across the nested functions (instead of directly accessing
    % object's properties) increases simulation speed.

    enum     = this.enum;     % Enumeration constants
    N_dim    = this.N_dim;    % Number of space dimensions
    N_p      = this.N_p;      % Number of particles
    T_emit   = this.T_emit;   % Emitter period between two particles, s (may be Inf)
    g        = [];            % Standard acceleration of Gravity (intensity), m/s^2
    gamma    = this.gamma;    %#ok<NASGU> Type of flow, 1 = laminar, 2 = turbulent
    imodel   = this.imodel;   % Interaction model
    k_s      = this.k_s;      % Spring coefficient
    k_d      = this.k_d;      % Spring damping factor
    k_e      = this.k_e;      % Impulse collision restitution (0 == max dissipation)
    k_p      = this.k_p;      % Position projection fraction (1 == max projection)
    mu_part  = this.mu_part;  %#ok<NASGU> Friction coefficient, particle-particle contact
    mu_sfc   = this.mu_sfc;   %#ok<NASGU> Friction coefficient, in contact with surface
    stepper  = this.stepper;  % Integrator method
    h        = NaN;           % Integrator time-step, s
    NT       = 0;             % Number of time-steps
    t_0      = this.t_0;      % Initial time, s
    t_f      = this.t_f;      % Final (simulation) time, s
    t        = this.t_0;      % Current simulation time, s
    t_real   = 0;             % Elapsed real time, s
    t_offset = 0;             % Difference between (toc) and (t), s
    t_idle   = 0;             % Idle (i.e. while sleeping) time, s
    t_incalc = 0;             % Time spent while solving equations, s
    t_indraw = 0;             % Time spent while rendering animation, s
    t_inupd  = 0;             % Time spent while updating plot data, s

    ZeroVec  = zeros( 1, N_dim );   % A helpfull 'zero' spatial vector
    
    % Particle state variables and properties
    M        = [];          % Mass
    R        = [];          % Radius; dead particles have NaN radius
    X        = [];          % Position
    V        = [];          % Velocity
    W        = [];          % Angular velocity
    FaceRgba = [];          % Face color (RGB + alpha)
    k_D      = [];          % Viscous drag coefficient
    k_M      = [];          % Magnus force coefficient
    k_W      = [];          % Angular velocity damping coefficient

    L        = [];          % System dimensions 3D: X,Y,Z, 2D: X,Z, 1D: Z
    LTick    = this.LTick;  % Characteristic length of the system

    XSFC     = [];          % System boundary surfaces' definint point
    NSFC     = [];          % System boundary surfaces' unit normal vector

    HF_z     = this.HF_z;   % Height field; in ndgrid format (x = rows, y = columns)
    HF_ends  = [];          % Height field x-y limits matrix [ xMin, xMax; yMin, yMax ]

    X_sink   = this.X_sink; % Sink position (only x, y coordinates; z is ignored)
    R_sink   = this.R_sink; % Sink radius
    
    track    = this.track;  % Structure with the trace of the system variables
    trace_T  = false;       % Flag indicating if time is traced
    trace_XP = false;       % Flag indicating if X, P, W and collision stats are traced
    trace_E  = false;       % Flag indicating if energy levels are traced

    % Animation related variables

    OH       = [];          % Graphic object handles
    trajPlot = [];          % Trajectory handles
    US_x     = [];          % Unit surface/circle
    US_y     = [];          % Unit surface/circle
    US_z     = [];          % Unit surface/circle
    pAsDots = false;        % Plot particles as dots if spheres are too small to be drawn

    mainFig  = this.mainFig;    % Main animation figure handle
    mainAxes = this.mainAxes;   % Main axes
    textInfo = this.textInfo;   % Text info handle
    miCancel = [];              % Cancel menu item handle
    movieObj = [];              % Movie object for capturing animations
    t_frame  = 1 / this.fps;    % Time per frame = reciprocal of the frames per second.
                                % Note that t_frame = Inf for fps = 0, which is OK.

	N_p_drawn  = 0;  % Number of particles that are painted last time
    t_r_redraw = 0;  % Scheduled real time for next animation frame, s
    t_e_redraw = 0;  % Scheduled simulation (elapsed) time for next animation frame, s
    wigglec    = 0;  % Wiggle counter used for wiggle stereoscopy

    v_fluid = this.v_fluid;      % Fluid velocity (e.g wind), m/s
    v_fluid(N_dim) = 0;          % z-component should be always 0
    v_fluid(N_dim+1:end) = [];   % Ensure correct number of dimensions
    this.v_fluid = v_fluid;      % Save back fixed value
    
    %=====================================================================================
    
    %%% Initialize simulation ------------------------------------------------------------

    initializeVariables ();
    initializeAnimation ();
    
    %%% Run simulation -------------------------------------------------------------------

    if ~this.Running
        return
    end

    mainLoop ();
    
    %%% Post-processing ------------------------------------------------------------------

    post_UpdateCaption ();
    post_CloseAnimation ();

    post_SaveData ();
    post_PlotEnergy ();

    %%% Return the results to the caller -------------------------------------------------

    if nargout
        varargout{1} = track;
    elseif ~this.Silent
        disp( 'Results:' ), disp( this  );
        disp( 'Track:'   ), disp( track );
    end

    %%% Indicate that the simulation is done ---------------------------------------------
    this.StopSimulation ();
    return

%=========================================================================================
%% NESTED FUNCTIONS ----------------------------------------------------------------------

%=========================================================================================
%% Initialization
% Create the data structures for the simulation and assigning them initial data.

function initializeVariables ()

    if ~this.Running
        return
    end

    %%% Assert the number of space dimensions --------------------------------------------

    N_dim = round( N_dim );
    assert( 1 <= N_dim && N_dim <= 3, ...
        '\n%s: The number of space dimensions must be either 1, 2 or 3', mfilename );

    % Setup acceleration of gravity as a vector and check number of dimensions
    % Assume gravity along the vertical axis (always denoted as Z-axis).

    switch N_dim
        case 1,  g = -this.g_n;  %#ok<SETNU>
        case 2,  g = [ 0, -this.g_n ];
        case 3,  g = [ 0, 0, -this.g_n ];
    end

    % Initialize 'L' as spatial vector holding dimensions of the system (box)
    % NOTE 1) The height is always the last dimension in L, i.e. L(end).
    % NOTE 2) Box always starts at [0,0,0] and ends at L
    % NOTE 3) If L(i) is negative, it is considered as 'non-limiting' box side
    %         i.e. the box side without the surface that collides with particles

    switch N_dim
        case 1,  L = [ this.BoxHeight ];
        case 2,  L = [ this.BoxWidth, this.BoxHeight ];
        case 3,  L = [ this.BoxWidth, this.BoxDepth, this.BoxHeight ];
    end

    %%% Initialize surfaces --------------------------------------------------------------

    % A surface is defined by some point on the surface XSFC, and 
    % the unit vector normal to the surface NSFC. Note that norm(NSFC) must be 1.
    % The following surfaces enclose the system in a box given by sizes in 'L'.
    % Note that negative values in L will not be considered as 'limiting' surfaces.

    XSFC = [ zeros( N_dim ); diag( L ) ];    % Points on boundary surfaces
    NSFC = [ eye( N_dim ); -eye( N_dim ) ];  % % Unit vectors on boundary surfaces

    % Remove surfaces that are negative, i.e. not limiting the box
    [infRows,~] = find( XSFC < 0 );
    infRows = unique( infRows );
    XSFC(infRows,:) = [];
    NSFC(infRows,:) = [];
    
    L = abs( L ); % Finally, make all dimensions positive

    %%% Initial state variables: particle mass, radius, position, velocity ---------------

    % Don't allow to add interactively particles in 3-D
    if this.ManualNp && N_dim == 3
        this.StopSimulation ();
        return
    end

    % Reset number of particles and traces, if particles should be placed interactivelly
    if this.ManualNp
        track = struct( 'n', 0 );
        N_p = 0;
    end

    % ---------------------------- Initialize particle status flag
    if ~isempty( this.M ) && all( [ N_p, N_dim ] <= size( this.M ) )
        % Use provided particle masses (reduced to our dimensions, if larger)
        M = WoP.diminishDimensions( this.M, N_p, N_dim );
    else
        % Use default particle masses
        M = this.m_p * ones( N_p, N_dim );
    end

    % ---------------------------- Initialize particle masses
    if ~isempty( this.M ) && all( [ N_p, N_dim ] <= size( this.M ) )
        % Use provided particle masses (reduced to our dimensions, if larger)
        M = WoP.diminishDimensions( this.M, N_p, N_dim );
    else
        % Use default particle masses
        M = this.m_p * ones( N_p, N_dim );
    end

    % ---------------------------- Initialize particle radii
    if ~isempty( this.R ) && N_p <= size( this.R, 1 )
        % Use provided particle radii (reduced to our dimensions, if larger)
        R = this.R(1:N_p,:);
    else
        % Use default particle radii
        R = this.r_p * ones( N_p, 1 );
    end

    % ---------------------------- Initial positions
    if ~isempty( this.X ) && all( [ N_p, N_dim ] <= size( this.X ) )
        % Use provided particle positions (reduced to our dimensions, if larger)
        X = WoP.diminishDimensions( this.X, N_p, N_dim );
    else
        X = zeros( N_p, N_dim ); % Allocate space

        %%% Pile up particles in a box (in a grid) with a small distances inbetween

        if N_p
            d_max = 2 * max( max( R ) ) + this.dx_emit; % Raster of the grid
            i_bound = floor( L / d_max ) - 1; % Number of particles per grid dimension
            i_np = prod( i_bound );      % Maximal number of particles in the grid
            i_cur = zeros( 1, N_dim );   % Current position in the grid, starting from 0

            % Set grid reference position (bottom-left = position of the first particle)
            if this.dx_emit > 0
                X(1,:) = this.dx_emit + ( 1 + rand( 1, N_dim ) ) * d_max/2;
            else
                X(1,:) = R(1) .* ones( 1, N_dim ); % put on the ground, if dx_emit == 0
            end
            % however, if a single particle, set it it's hight is fully random
            if N_p == 1
                X(end,:) = ( 0.5 + 0.5 * rand(1) ) * ( L(end) - d_max );
            end
        end

        % Place particles snapped to the grid
        for i = 1:N_p
            X(i,:) = X(1,:) + i_cur * d_max;

            i_cur(end) = i_cur(end) + 1; % Increment position in the grid along z-axis
            if i >= i_np
                % Grid is full; continue filling along z-axis...
                % (this also suppresses filling along non-existing dimensions)
                continue
            elseif i_cur(end) >= i_bound(end)
                % Reached z-axis bound; now continue along x-axis...
                i_cur(end) = 0;
                i_cur(1) = i_cur(1) + 1;
                if i_cur(1) >= i_bound(1)
                    % Reached x-axis bound; now continue on y-axis...
                    i_cur(1) = 0;
                    i_cur(2) = i_cur(2) + 1;
                end
            end
        end
    end

    % ---------------------------- Initialize particle face colors

    if ~isempty( this.FaceRgba ) && N_p <= size( this.FaceRgba, 1 )
        % Use provided face colors (reduced to our N_p, if larger)
        FaceRgba = this.FaceRgba(1:N_p,:);
    else
        % Generate particle face colors: pastel, with particle# as hue
        % Note that 4-th component of FaceRgba is alpha (transparency)
        FaceRgba = zeros( N_p, 4 );
        for i = 1:N_p
            FaceRgba(i,:) = this.randomFaceColor( ( i - 1 ) / N_p );
        end
    end

    % ---------------------------- Initial velocities
    if ~isempty( this.V ) && all( [ N_p, N_dim ] <= size( this.V ) )
        % Use provided particle velocities (reduced to our dimensions, if larger)
        V = WoP.diminishDimensions( this.V, N_p, N_dim );
    else
        V = zeros( N_p, N_dim ); % Allocate space

        if this.g_n == 0
            % Generate random velocities without, if no gravity
            V = bsxfun( @times, -1 + 2 * rand( N_p, N_dim ), L );
        elseif N_dim >= 2
            % Induce very small velocity disturbances in horizontal plane
            V(:,1:end-1) = -this.dv_emit + 2 * this.dv_emit * rand( N_p, N_dim-1 );
        end
    end

    % ---------------------------- Initial angular velocities
    if ~isempty( this.W ) && all( [ N_p, 3 ] <= size( this.W ) )
        % Note that angular velocity is always 3-D spatial vector
        W = WoP.diminishDimensions( this.W, N_p, 3 );
    else
        W = zeros( N_p, 3 ); % Allocate space
    end

    %%% Initialize calculated quantities -------------------------------------------------

    rho = this.rho;  % Cache some values to be more verbose
    C_D = this.C_D;
    C_M_at_W = this.C_M / this.W_M;  % See note for C_M/W_M in WoP.m
    C_W = this.C_W;

    A   = pi * R .^ 2;                   % Cross section area
    k_D = 0.5 * rho * A * C_D;           % Drag coefficent, kg/m
    k_M = 0.5 * rho * A * C_M_at_W;      % Magnus force coefficient, kg.s/m
    k_W = 0.5 * rho * A * C_W ./ M(:,1); % Magnus force coefficient, kg.s/m

    %%% Setup the integrator variables ---------------------------------------------------

    assert( stepper == enum.ForwardEuler      ...
         || stepper == enum.SemiImplicitEuler ...
         || stepper == enum.Leapfrog, ...
        '\nUnknown integration method; stepper = %d', stepper );

    if stepper == enum.ForwardEuler && imodel == enum.ImpulseModel
        throw( MException( 'WoP:Simulation:feim', ...
            'Forward Euler integrator is not supported for the Impulse Model' ...
        ) );
    end

    % Time-step depends on collision model:
    switch imodel
        case enum.SpringModel,  h = this.h_soft;  % Time-step size for linear spring model
        case enum.ImpulseModel, h = this.h_hard;  % Time-step size for impulse model
    end
    assert( h > 0 );
    this.h = h;  % Remember selection

    % Number of time steps must be integer
    NT = ceil( ( t_f - t_0 ) / h );

    %%% Compute derived quantities -------------------------------------------------------

    if this.ManualNp
        return
    end

    result = WoP_Solver;

    %%% Allocate vectors/matrices for traced variables -----------------------------------

    % Be ware that tracing from one simulation can be continued into the next!

    % Do we need to track particle trajectory?
    % Yes, if the user wants to or we are already tracing ('X' field exists).
    trace_XP = this.TraceVars || nargout >= 1 || isfield( track, 'X' ) ...
        || this.PlotPath || this.PlotContacts || this.PlotMomentum || this.PlotHeight ...
        || this.PlotAngularW;

    % Do we need to track energy of the system over time?
    % Yes, if the user wants to or we are already tracing ('E' field exists).
    trace_E  = this.TraceVars || isfield( track, 'E' ) || this.PlotEnergy;
    
    % Do we need to keep track of elapsed time?
    trace_T = trace_XP || trace_E;

    % Why do we need to track time? Coludn't it be generated as t = t_0 + n * h?
    % Yes, but only for the current simulation, as the time-step h may be different 
    % accross several simulations!

    % In case that earlier track was for less particles, fill
    % the missing variable traces with default values
    if trace_XP && isfield( track, 'X' )
        oldnp = size( this.track.X, 2 );
        if oldnp < N_p
            this.track.X( :, oldnp + 1 : N_p, : ) = NaN;
            this.track.P( :, oldnp + 1 : N_p, : ) = 0;
            this.track.W( :, oldnp + 1 : N_p, : ) = 0;
        end
    end

    % Keep track of the current time-step number (where track.n = 1 corresponds to t_0)
    if track.n
        % Proceed with n from the last simulation
    else
        track.n = 1; % Initialize time-step number
    end

    % Alocate/expand vector holding elapsed time
    if trace_T
        if ~isfield( track, 'T' )
            track.T = NaN * ones( track.n + NT, 1 );
        else
            track.T = cat( 1, track.T, NaN * ones( NT, 1 ) );
        end
        % Intial value
        track.T( track.n ) = t_0;
    end

    % Alocate/expand matrices holding position and linear momentum trajectories
    if trace_XP
        if ~isfield( track, 'X' )
            % Allocate
            track.X    = NaN * ones( [ track.n + NT, N_p, N_dim ] );
            track.P    = zeros( [ track.n + NT, N_p, N_dim ] );
            track.W    = zeros( [ track.n + NT, N_p, 3 ] );
            track.colc = zeros( track.n + NT, 3 );
        else
            % Expand
            track.X    = cat( 1, track.X,    NaN * ones( [ NT, N_p, N_dim ] ) );
            track.P    = cat( 1, track.P,    zeros( [ NT, N_p, N_dim ] ) );
            track.W    = cat( 1, track.W,    zeros( [ NT, N_p, 3 ] ) );
            track.colc = cat( 1, track.colc, zeros( NT, 3 ) );
        end
        % Intial values
        track.X( track.n, :, : ) = X;
        track.P( track.n, :, : ) = M .* V;
        track.W( track.n, :, : ) = W;
    end

    % Allocate/expand solutionvectors holding energy trace and track initial values
    if trace_E
        if ~isfield( track, 'E' )
            track.E = NaN * ones( track.n + NT, 3 ); 
        else
            track.E = cat( 1, track.E, NaN * ones( NT, 3 ) );
        end
        % Intial energy values are calculated by calling solver without arguments
        track.E( track.n, : ) = result.E;
    end

end % initializeVariables

%=========================================================================================
%% Initialization of the Animation, if enabled
% Initialize the main figure that will be used for the animation, together with graphics 
% objects representing particles and their trajectories.

function initializeAnimation ()

    % Proceede only if still running and animation is enabled or manual placing of 
    % particles is enabled
    if ~( this.Running && ( this.Animate || this.ManualNp ) )
        return
    end

    if max( L ) / this.LTick > 20
        this.LTick = max(L)/20;
    end

    % Get the size of the screen
    set( 0, 'Units', 'Pixels' );
    screen.Rect = get( 0, 'ScreenSize' ); % gets [ left, bottom, width, height ]
    screen.L = screen.Rect(1);  screen.B = screen.Rect(2);
    screen.W = screen.Rect(3);  screen.H = screen.Rect(4);

    % Calculate paper dimensions for the figure (height and width) which
    % depends whether the picture is flat or 3-D
    switch N_dim
        case 1,  screen.box_W = 2 * LTick;
                 screen.box_H = L(1);
        case 2,  screen.box_W = L(1);
                 screen.box_H = L(2);
        case 3,  screen.box_W = ( L(1) + L(2) ) * 0.86;
                 screen.box_H = L(3) + max( L(1), L(2) ) * 0.50;
    end

    % The list of axes properties retained between simulations.
    % (Don't forget to put 'prop*Mode' after 'prop*' when maintaining the list!)

    retainedView = [];
    retainedAxesPropsValues = {};  % If empty, props will be initialized to defaults

    retainedAxesProps = { ...
        'Box';              'Projection';
        'DataAspectRatio';  'DataAspectRatioMode'; 
        'CameraPosition';   'CameraPositionMode'; 
        'CameraTarget';     'CameraTargetMode'; 
        'CameraUpVector';   'CameraUpVectorMode'; 
        'CameraViewAngle';  'CameraViewAngleMode';
        'TickLength';     'TickDir';        'TickDirMode';
        'XColor';         'YColor';         'ZColor';
        'XDir';           'YDir';           'ZDir';
        'XGrid';          'YGrid';          'ZGrid';
        'XLim';           'YLim';           'ZLim';
        'XLimMode';       'YLimMode';       'ZLimMode';
        'XTick';          'YTick';          'ZTick';
        'XTickMode';      'YTickMode';      'ZTickMode';
        'XTickLabel';     'YTickLabel';     'ZTickLabel';
        'XTickLabelMode'; 'YTickLabelMode'; 'ZTickLabelMode';
    };

    % Open an existing or create new main figure
    
    if ~isempty( mainAxes ) && ishandle( mainAxes )

        % If we got only axes handle, set mainFig as figure holding the axes
        mainFig = get( mainAxes, 'Parent' );
        if ~strcmp( get( mainFig, 'Type' ), 'figure' )
            mainFig = get( mainFig, 'Parent' );
        end
        if ~strcmp( get( mainFig, 'Type' ), 'figure' )
            error( 'Parent or grand-parent of the mainAxes must be a figure.' );
        end

        set( 0, 'CurrentFigure', mainFig );

        if ~isequal( get( mainAxes, 'UserData' ), 'WoP' )
            % Not our axes. Just clear axes from objects.
            cla( mainAxes );
        else
            % Our axes. Remember axes and view properties
            retainedAxesPropsValues = get( mainAxes, retainedAxesProps );
            [retainedView(1),retainedView(2)] = view( mainAxes );

            % Clear axes from objects
            cla( mainAxes );

            % Restore axes and view properties
            view( mainAxes, retainedView )
            set( mainAxes, retainedAxesProps, retainedAxesPropsValues );
        end

    elseif ~isempty( mainFig ) && ishandle( mainFig )

        % Erase previous contents of the figure
        set( 0, 'CurrentFigure', mainFig );
        mainAxes = gca;

        % Remember axes and view properties
        retainedAxesPropsValues = get( mainAxes, retainedAxesProps );
        [retainedView(1),retainedView(2)] = view( mainAxes );

        % Clear figure from objects (including axes)
        clf( mainFig );

        % Restore axes and view properties
        mainAxes = gca; % Creates new axes
        view( mainAxes, retainedView )
        set( mainAxes, retainedAxesProps, retainedAxesPropsValues );
    else
        % Create a new main figure
        this.mainFig = figure( ...
            'NumberTitle', 'off', ...
            'Name', [ this.Title, ' (', this.ID, ')' ], ... % window title
            'FileName', [ 'lab3_fig1_', this.ID, '.fig' ], ... % name of the FIG-file 
            'Color', [ 0.97, 0.97, 0.9 ] ...
        );

        % Normalize paper height/width to max. allowed figure width (in centimeters)
        screen.box_maxWH = max( screen.box_W, screen.box_H );
        screen.pap_W = 3 + ( this.MaxFigSize - 3 ) * screen.box_W / screen.box_maxWH;
        screen.pap_H = 3 + ( this.MaxFigSize - 3 ) * screen.box_H / screen.box_maxWH;
        screen.pap_maxWH = max( screen.pap_W, screen.pap_H );

        % Determine windows borders and reduce screen size
        screen.borders = get( this.mainFig, 'OuterPosition' ) ...
                       - get( this.mainFig, 'Position' );

        % Reduce screen size for windows borders and windows taskbar height
        screen.W = screen.W - screen.borders(3);
        screen.H = screen.H - screen.borders(4) - 40;
        screen.addLeft = 1  - screen.borders(1);
        screen.addBot  = 40 - screen.borders(2);

        % Printer to screen conversion factor (pixels per centimeter), so the whole 
        % printer figure fits in PctFullScreen of the screen
        screen.dpcm = min( screen.H / screen.pap_H, screen.W / screen.pap_W ) ...
                    * this.PctFullScreen;

        % Configure screen size and position (position is centered with slight offset)
        screen.fig_W = screen.pap_W * screen.dpcm; % Figure width on screen
        screen.fig_H = screen.pap_H * screen.dpcm; % Figure height on screen
        screen.fig_L = screen.addLeft + ( screen.W - screen.fig_W ) * 0.3; % left pos.
        screen.fig_B = screen.addBot  + ( screen.H - screen.fig_H ) * 0.7; % bottom pos.

        set( this.mainFig, ...
            'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
            'PaperSize', [ screen.pap_W, screen.pap_H ], ...
            'PaperPosition', [ 0, 0, screen.pap_W, screen.pap_H ], ...
            'Position', [ screen.fig_L, screen.fig_B, screen.fig_W, screen.fig_H ] ...
        );

        mainFig  = this.mainFig;

        % Configure toolbar visibility
        if N_dim == 3
            set( mainFig, 'Toolbar', 'none' ) % hide standard toolbar
            cameratoolbar( mainFig, 'Show' )  % show camera toolbar
        end

        % Create a new axes
        mainAxes = gca;
    end

    % Set flicker-free rendering
    set( mainFig, 'DoubleBuffer', 'on' );

    % Add 'Stop Simulation' to the menu, if not manually placing particles or in Golf
    if ~this.ManualNp && ~this.GolfTweaks
        set( mainFig, 'UserData', this );
        miCancel = uimenu( mainFig, 'Label', '&Stop', ...
            'Callback', 'StopSimulation( get( gcbf, ''UserData'' ) );', ...
            'ForegroundColor', [ 0.7, 0, 0 ] ...
        );
    end

    % Add context menu items to export figure as EPS/EMF/PNG
    cxMenu = uicontextmenu( 'Parent', mainFig );
    uimenu( cxMenu, 'Label', 'Export as EPS', 'Callback', [ ...
        'WoP.ExportFigure( gcbf, ''-depsc2'', ', ...
        '''lab3_fig1_', this.ID, ''' );' ...
    ]);
    uimenu( cxMenu, 'Label', 'Export as PNG', 'Callback', [ ...
        'WoP.ExportFigure( gcbf, ''-dpng'', ', ...
        '''lab3_fig1_', this.ID, ''' );' ...
    ]);
    uimenu( cxMenu, 'Label', 'Export as EMF', 'Callback', [ ...
        'WoP.ExportFigure( gcbf, ''-dmeta'', ', ...
        '''lab3_fig1_', this.ID, ''' );' ...
    ]);
    % Add context menu item to export annotation as EPS
    uimenu( cxMenu, 'Label', 'Export Annotation', 'Separator', 'on', ...
        'Callback', [ ...
            'WoP.ExportAnnotation( gcbf, ', ...
            '''figCaption'', ''-depsc2'', ''', ...
            [ 'lab3_info_', this.ID ], ''' );' ...
        ] ...
    );
    set( mainFig, 'UIContextMenu', cxMenu );
    set( mainAxes, 'UIContextMenu', cxMenu );

    % Initialize axes and view (if not already initialized) ------------------------------
    
    set( mainAxes, 'UserData', 'WoP' );  % Mark axes as ours

    if isempty( retainedView )
        view( 146, 15 )
    end

    if isempty( retainedAxesPropsValues )
        % Default values for a new figure
        axis equal; % Sets equal scaling on all axes
        tck = LTick;   % Distance between ticks
        gridColor = [ 0.5, 0.8, 0.8 ];  % dark cian

        switch N_dim   % Axes, N_dim dependent part:
          case 1
            % Axes limits
            axis( [ -tck, tck, -0.2, L(1) ] );
            % Ticks
            set( mainAxes, 'XTick', -tck:tck:tck, 'YTick', 0:tck:L(1) );
            % Grid
            grid on;
            % Projection
            set( mainAxes, 'Projection', 'orthographic' );
            view( 2 );
          case 2
            % Axes limits
            axis( [ 0, L(1), -0.2, L(2) ] );
            % Ticks
            set( mainAxes, 'XTick', 0:tck:L(1), 'YTick', 0:tck:L(2) );
            % Grid
            grid on;
            % Projection
            set( mainAxes, 'Projection', 'orthographic' );
            view( 2 );
          case 3
            % Axes
            axis( [ 0, L(1), 0, L(2), 0, L(3) ] );
            set( mainAxes, 'XTickLabel', '', 'YTickLabel', '', 'ZTickLabel', '' );
            % Ticks
            set( mainAxes, 'XTick', 0:tck:L(1), 'YTick', 0:tck:L(2), 'ZTick', 0:tck:L(3) );
            set( mainAxes, 'TickDir', 'out' );  % Ticks outside
            % Grid
            set( mainAxes, 'XColor', gridColor, 'YColor', gridColor, 'ZColor', gridColor );
            grid on;
            box on;
            % Projection
            set( mainAxes, 'Projection', 'perspective' );
            if this.GolfTweaks
                grid off;
                set( mainAxes, ...
                    'TickLength', [ 0, 0 ], ...
                    'CameraPositionMode',   'manual', ...
                    'CameraPosition',       [ 0.5 * L(1), -1 * L(2), 1.4* L(3) ], ...
                    'CameraTargetMode',     'manual', ...
                    'CameraTarget',         [ 0.5 * L(1), 0.5 * L(2), 0.3 * L(3) ], ...
                    'CameraUpVectorMode',   'manual', ...
                    'CameraUpVector',       [ 0, 0, 1 ], ...
                    'CameraViewAngleMode',  'manual', ...
                    'CameraViewAngle',      20 ...
                );
            end
        end
    end

    % (Re)define axes labels

    switch N_dim
        case 1
            ylabel( '{\itz} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
        case 2
            xlabel( '{\itx} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
            ylabel( '{\itz} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
    end

    % Plot the contents of the figure ----------------------------------------------------

    hold on; % Retain subsequent plots in figure

    % Init paticle visual representation template: either a sphere or a circle

    if N_dim == 3   % Render particles as spheres always when plotting in 3D
        this.PAsSpheres = true;
    end

    % If particle size is too small turn-on showing particles as dots
    if 800 * this.r_p / max( L ) < 1
        % pAsDots = true;
    end

    if pAsDots
        US_x = 0;
        US_y = 0;
        US_z = 0;
    elseif this.PAsSpheres
        % Generates surface points for a unit sphere
        if N_dim == 3
            [ US_x, US_y, US_z ] = sphere( this.N_usf );
        else
            % Swapped Y/Z axies when plotting in 1- or 2-D
            [ US_x, US_z, US_y ] = sphere( this.N_usf );
        end
        % Also, turn-on OpenGL if ploting surfaces
        set( mainFig, 'Renderer', 'OpenGL' );
    else
        % Generates edge points for a unit circle
        ang = linspace( 0, 2 * pi, this.N_uck );
        [ US_x, US_y, US_z ] = deal( cos(ang), sin(ang), zeros(size(ang)) );
    end

    % Use the light source when displaying particles as spheres

    if this.PAsSpheres || N_dim == 3
        if this.GolfTweaks
            light( 'Position', [ 10, -15, 5 ], 'Style', 'infinite' );
        else
            light( 'Position', [ -10, 10, 4 ], 'Style', 'infinite' );
        end
    end

    % Generate ground surface as height field. Creates several peaks and valleys 
    % using the gaussian exponential function accumulating these to the height field.

    switch N_dim
        case 1, HF_ends = [];
        case 2, HF_ends = [ 0, L(1); 0, 0 ];
        case 3, HF_ends = [ 0, L(1); 0, L(2) ];
    end

    if ~this.EnableHF  
        HF_z = []; % Disable heightfield
    elseif isempty( HF_z )
        switch N_dim
            case 1
                HF_z = []; % No height fields in 1-D
            case 2
                HF_z = WoP.GenerateHeightField( HF_ends, 50, 0.5 * LTick );
            case 3
                HF_z = WoP.GenerateHeightField( HF_ends, 50, 1 * LTick );
        end
    elseif N_dim == 1  % turn-off height field in 1-D
        HF_z = [];
    end

    % Now we can draw ground level, either as a height field or common line/surface

    if ~isempty( HF_z ) 
        % Draw the height field either as a horizontal black thick line in 1-D
        % or a dark-green thick curve in 2-D or a half-transparent green surface in 3-D
        switch N_dim
          case 1
            plot( [ -LTick, LTick ], [ 0, 0 ], '-k', 'LineWidth', 1.5 );
          case 2
            hfx = linspace( HF_ends(1,1), HF_ends(1,2), size(HF_z,1) );
            plot( hfx, HF_z', '-', 'DisplayName', 'Ground Surface', ...
                'Color', [ 0, 0.5, 0 ], 'LineWidth', 1.5 );
          case 3
            [ hfx, hfy ] = meshgrid( ... % Recover meshgrid for the height field
                linspace( HF_ends(1,1), HF_ends(1,2), size(HF_z,1) ), ...
                linspace( HF_ends(2,1), HF_ends(2,2), size(HF_z,2) ) );
            surf( hfx, hfy, HF_z', ... % Tranposed height field surface (from ndgrid)
                'DisplayName', 'Ground Surface', ...
                'EdgeColor', [ 0, 1.0, 0 ], 'EdgeLighting', 'gouraud', ...
                'FaceColor', [ 0.2, 0.9, 0.2 ], 'FaceLighting', 'gouraud' ...
            );
        end
    else
        % Draw the ground level either as horizontal thick black line in 1- or 2-D
        % or a half-transparent yellowish plane crossing point [0,0,0] in 3-D
        switch N_dim
          case 1
            plot( [ -LTick, LTick ], [ 0, 0 ], '-k', 'LineWidth', 1.5 );
          case 2
            plot( [ -LTick, L(1) + LTick ], [ 0, 0 ], '-k', 'LineWidth', 1.5 );
          case 3
            fill( [ 0, L(1), L(1), 0, 0 ], [ 0, 0, L(2), L(2), 0 ], [ 1, 1, 0 ], ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none' ); % Half-transparent without edges
        end
    end

    % Remove bottom surface in case that we have a height field

    if ~isempty( HF_z ) && N_dim > 1
        XSFC(N_dim,:) = [];
        NSFC(N_dim,:) = [];  %#ok<SETNU>
    end

    % Save the height field into the object's store

    if this.EnableHF
        this.HF_z = HF_z;
    end

    % Draw sink as black rectangle in 2-D or black cylinder in 3-D
    
    if ~isempty( X_sink ) && ~isempty( R_sink )
        switch N_dim
          case 2
            rectangle( 'Position', [  X_sink(1) - R_sink, 0,  2 * R_sink, X_sink(2) ], ...
                'FaceColor', 'k', 'EdgeColor', 'k' );
          case 3
            [ c_x, c_y, c_z ] = cylinder;
            surf( ...
                c_x * R_sink + X_sink(1), c_y * R_sink + X_sink(2), c_z * X_sink(3), ...
                'EdgeLighting', 'gouraud', 'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud', 'FaceColor', 'k' ...
            );
            plot3( [ X_sink(1), X_sink(1) ], [ X_sink(2), X_sink(2) ], ...
                  [ 0, 1.5*X_sink(3) ], '-k', 'LineWidth', 1.5 ...
            );
        end
    end

    % Capture movie; prepare the new file and create an animation.

    if this.CaptureMovie
        movieObj = VideoWriter([ 'lab3_', this.ID, '.avi' ]);
        movieObj.FrameRate = this.fps;
        movieObj.Quality = 75;
        open( movieObj );
    end

    %%% Initialize graphics object handles associated to particles -------------------

    OH = zeros( N_p, 1 ); % Allocate space for particle graphics object handles

    % Allocate space for trajectory plot handles
    if this.PlotPath
        trajPlot = zeros( N_p, 1 ); 
    end

    % Finaly, create graphical objects representing particles
    for i = 1:N_p
        initParticle_Visual( i );
    end

    % Display velocity vector and particle state variables, if not manual placing
    % mode and if AnnotInit is set or initial and final times are equal and number
    % of particles is less or equal 8).
    if ~this.ManualNp && ( this.AnnotInit || ( t_0 == t_f && this.N_p <= 8 ) )
        for i = 1:N_p
            annotateParticle( i, struct( 'm', M(i,:), 'r', R(i,:), ...
                'x', X(i,:), 'v', V(i,:), 'w', W(i,:), 'color', FaceRgba(i,:) ) );
        end
    end

    % Display simulation status in the top region of the figure
    if isempty( textInfo ) || ~ishandle( textInfo )
        textInfo = text( 0.5, 1.01,  '', 'Tag', 'WoP_textInfo', ...
            'Units', 'normalized', ...
            'FontName', 'Helvetica', 'FontSize', 10, 'Color', [ 0, 0, 0.9 ], ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none', 'EdgeColor', 'none' ...
        );
    end

    % Popup figure
    figure( mainFig );

    %%% Interactivelly add new particles ---------------------------------------------

    if this.ManualNp

        addNewParticles_User ();

        this.ManualNp = false;  % we're done; turn off the 'manual mode' flag
        this.N_p = N_p;         % update number of particles

        if N_p == 0 % If no particles were added
            delete( mainFig );   % close the window
        else
            % Save the particle state variables ...
            this.M_0 = M;  % initial masses
            this.R_0 = R;  % initial radii
            this.X_0 = X;  % initial positions
            this.V_0 = V;  % initial velocities
            this.W_0 = W;  % initial velocities
            this.M   = M;  % final masses
            this.R   = R;  % final radii
            this.X   = X;  % final positions
            this.V   = V;  % final velocities
            this.W   = W;  % final velocities
            % ... and particle colors:
            this.FaceRgba_0 = FaceRgba;  % initial face colors
            this.FaceRgba   = FaceRgba;  % final face colors
        end

        this.StopSimulation ();
    end

    drawnow;  % Flush event queue and update figure window
    set( 0, 'CurrentFigure', mainFig );  % Go back to mainFig

end % initializeAnimation

%=========================================================================================
%% Sanity-check of the simulation parameters

function sanityCheck ()

    assert( isequal( size( M ), [ N_p, N_dim ] ) );
    assert( isequal( size( R ), [ N_p, 1 ] ) );
    assert( isequal( size( X ), [ N_p, N_dim ] ) );
    assert( isequal( size( V ), [ N_p, N_dim ] ) );
    assert( isequal( size( W ), [ N_p, 3 ] ) );
    assert( isequal( size( FaceRgba ), [ N_p, 4 ] ) );
    assert( all( M(:) > 0 ) );
    assert( all( R(:) >= 0 | isnan( R(:) ) ) );

end % function sanityCheck

%=========================================================================================
%% Main Simulation Loop

function mainLoop ()

    sanityCheck ();

    % If not animating particles, open the progress bar
    if ( ( ~this.Animate && t_f > t_0 ) || isinf( t_frame ) )
        t_waitbInc = ( t_f - t_0 ) / 20;   % Waitbar update each 5% of time
        t_waitbUpdate = t_0 + t_waitbInc;  % Schedule next update
        waitbh = waitbar( 0, sprintf( 't = %.3f s', t ), ...
            'Name', [ 'Running ', this.Title ], ...
            'WindowStyle', 'modal', 'Color', [ 0.925, 0.914, 0.847 ], ...
            'UserData', this, ...
            'CreateCancelBtn', [ 
                'try, StopSimulation( get( gcbf, ''UserData'' ) ); catch, end, ', ...
                'delete( gcbf )' ...
            ] ...
        );
    else
        % There is no progress bar = no need for updates
        t_waitbUpdate = Inf;
    end

    % Start measuring elapsed (real) time and remember current simulation time (t)
    % could be compared later with (t_offset + toc).
    tic;
    t_start    = 0;   % for elapsed (real) time statistics
    t_offset   = t;   % real time = t_ffset + toc
    
    % Initialize next time to draw animation frame
    t_e_redraw = t;   % scheduled next redraw, ref to simulated time
    t_r_redraw = t;   % scheduled next redraw, ref to real time

    if ~this.InRealTime
        t_r_redraw = Inf;  % Cancel real-time updates if not in real time
    end
    if isinf( t_frame )
        t_r_redraw = Inf;  % Cancel real-time updates if fps = 0
        t_e_redraw = Inf;  % Cancel elapsed-time updates if fps = 0
    end

    % Initialize next time to emit a new particle
    t_nextEmit = t_0 + T_emit - h;  

    % Remember the initial state variables
    if t_0 == 0
        this.M_0 = M;  % save initial masses and ...
        this.R_0 = R;  % initial radii
        this.X_0 = X;  % initial positions
        this.V_0 = V;  % initial velocities
        this.W_0 = W;  % initial angular velocities
        this.FaceRgba_0 = FaceRgba;  % initial face colors
    end

    solution.N_p = N_p;

    for n = 1 : NT

        %---------------------------------------------------------------------------------
        %% » MEX function WoP_Solver does:
        %  - Collision detection and response
        %  - Calculation of internal and external forces
        %  - Integration of ODE of motion
        %  - Computation of derived quantities
        %---------------------------------------------------------------------------------

        solution = WoP_Solver( this.IgnorePP, this.KillInactive, this.KillOutOfBox );

        %---------------------------------------------------------------------------------
        % Solver returns a structure with fields: t, X, V, W, R, E, P, ccStat

        t = solution.t;
        X = solution.X;
        V = solution.V;
        W = solution.W;
        R = solution.R;

        %---------------------------------------------------------------------------------
        %% » Trace selected variables like position and momentum trajectory, or energies

        % Update the current track-record number, where track.n = 1 corresponds to t_0
        % i.e. track.n = n + 1 corresponds to t_n
        track.n = track.n + 1;

        % Save variables
        if trace_T
            track.T   ( track.n       ) = t; 
        end
        if trace_XP
            track.X   ( track.n, :, : ) = solution.X;
            track.P   ( track.n, :, : ) = solution.P;
            track.W   ( track.n, :, : ) = solution.W;
            track.colc( track.n, :    ) = solution.ccStat;
        end
        if trace_E
            track.E   ( track.n, :    ) = solution.E;
        end

        % Accumulate time spent in solving eqs.
        t_stop = toc;
        t_incalc = t_incalc + t_stop - t_start;
        t_start = t_stop;

        %---------------------------------------------------------------------------------
        %% Update progress bar, if scheduled

        if t >= t_waitbUpdate
            t_waitbUpdate = t + t_waitbInc;  % Schedule next update
            waitbar( ( t - t_0 ) / ( t_f - t_0 ), waitbh, sprintf( 't = %.3f s', t ) );
        end

        %---------------------------------------------------------------------------------
        %% Render graphis in real time

        if this.Animate
            % Determine what's a real time
            if this.CaptureMovie
                t_real = t;   % real time = simulated time
            else
                t_real = t_offset + toc; % real time = elapsed real time + initial offset
            end
            % Update particles' and trajectories' positions: either if scheduled, 
            % or if N_p has changed since last redraw
            if  t >= t_e_redraw || t_real >= t_r_redraw || N_p ~= N_p_drawn
                renderGraphics ();
            end
        end

        %---------------------------------------------------------------------------------
        %% Remove inactive particles, if enabled

        if solution.N_p ~= N_p && this.RemoveKilled
            removeParticle( find( isnan(R) ) );
        end

        %---------------------------------------------------------------------------------
        %% Stop the simulation, if there are no active particles or if stopped

        if ~this.Running || solution.N_p <= 0
            break;
        end

        %---------------------------------------------------------------------------------
        %% Emit a new particle, if scheduled

        if t >= t_nextEmit
            t_nextEmit = t_nextEmit + T_emit - h; % Setup next time to emit a new particle
            emitNewParticle_HiZ( track.n );
        end

        % Accumulate time spent in drawing
        t_stop = toc;
        t_indraw = t_indraw + t_stop - t_start; 
        t_start = t_stop;

    end % of the Main Simulation Loop ----------------------------------------------------

    this.Running = false;
    t_real = toc;  % Get elapsed real time

    % Save the resulting number of active particles
    this.N_active = solution.N_p; 

    %%% Render the last frame and close the progress bar (if exists)

    if NT ~= 0 && this.Animate
        renderGraphics ();
    end

    if exist( 'waitbh', 'var' ) && ishandle( waitbh )
        delete( waitbh );
    end

    %% Show the reason why we are here...

    if solution.N_p <= 0
        reason = 'STOPPED';     % There were no more active particles
    elseif t >= t_f
        reason = 'COMPLETED';   % The simulation has been successfully completed
    else
        reason = 'CANCELLED';   % The simulation has been cancelled
    end

    fprintf( '\n>>> %s at n = %d, t = %g s, t_f = %g s, NT = %d\n\n', ...
        reason, track.n - 1, t, t_f, NT );

end % function mainLoop

%-------------------------------------------------------------------------------------
%% » Real-Time Graphics Rendering

function renderGraphics ()

    N_p_drawn = N_p; % Update number of drawn particles

    % If we are here because redrawing was scheduled earlier
    if  t >= t_e_redraw || t_real >= t_r_redraw

        % Schedule next updates
        t_e_redraw = t_e_redraw + t_frame;
        t_r_redraw = t_r_redraw + t_frame;

        % Calculate whether to pause if rendering is too fast or not.
        % Time to sleep = simulation elapsed time - real elapsed time
        t_pause = t - t_real;

        % If too fast, pause so the simulation is in real-time
        % If too slow, increase frame period porportionally to delay but max 100 ms
        if this.InRealTime && ~isinf( t_frame )
            if t_pause > 0
                t_idle = t_idle + t_pause; % Update idle time statistics
                pause( t_pause - 0.005 );   % Reduce 5 ms for drawnow
                % fprintf( '%8.3f %8.3f  Late  %6.3f\n', t, t_real, t_pause );
            elseif t_pause < 0
                % fprintf( '%8.3f %8.3f  Late  %6.3f\n', t, t_real, t_pause );
                t_e_redraw = t + t_frame;
                t_r_redraw = t_r_redraw + min( -t_pause, 0.1 );
            end
        end
    end

    % Set main figure as current (without figure pop-up)
    set( 0, 'CurrentFigure', mainFig );

    t_start = toc;

    for i = 1:N_p

        % Skip particles without handles
        if ishandle( OH(i) )
            r = R(i);
            if isnan(r)
                % Remove inactive particles (with NaN radius)
                if ~this.GolfTweaks
                    delete( OH(i) );
                end
            else
                % Update particle postion, either as 3D-surface or 2D-polygon
                switch N_dim
                    case 1,  c = [ 0, X(i), 0 ];
                    case 2,  c = [ X(i,:), 0 ];
                    case 3,  c = X(i,:);
                end
                set( OH(i), 'XData', r * US_x + c(1), ...
                            'YData', r * US_y + c(2), ...
                            'ZData', r * US_z + c(3) );
            end
        end

        % Update particle trajectory, either as 2- or 3-D curve.
        % (Note that it's ridicilous to plot trajectories in 1-D.)
        if this.PlotPath
            switch N_dim
                case 2
                    set( trajPlot(i), ...
                        'XData', track.X( 1:track.n, i, 1 ), ...
                        'YData', track.X( 1:track.n, i, 2 ) );
                case 3
                    set( trajPlot(i), ...
                        'XData', track.X( 1:track.n, i, 1 ), ...
                        'YData', track.X( 1:track.n, i, 2 ), ...
                        'ZData', track.X( 1:track.n, i, 3 ) );
            end
        end
    end

    % Update status info

    set( textInfo, 'String', sprintf( 'N_p = %d,  t = %.3f s', N_p, t ) );

    t_inupd = t_inupd + toc - t_start; % Accumulate time spent while updating plot data

    %%% Wiggle stereoscopy, only if plotting in 3-D and if it's enabled

    if N_dim == 3 && this.Wiggle3D
        dtheta = 0.1; % smooth orbit camera
        % Wiggle camera left/right
        wigglec = mod( wigglec + 1, 20 );
        if wigglec >= 10
            camorbit( dtheta, 0 ); % half-time orbit camera right
        else
            camorbit( -dtheta, 0 ); % half-time orbit camera left
        end
    end

    %%% Screen update

    if this.CaptureMovie
        writeVideo( movieObj, getframe( mainAxes ) );
    else
        drawnow;
    end

end % function renderGraphics

%% Create common caption/annotation for the figures (in TeX format) 

function post_UpdateCaption ()

    % Prepare info about interaction model
    switch imodel
      case this.enum.ImpulseModel
        imodel_info = [ this.verb.imodel{imodel}, ...
          ': {\ite} = ',   sprintf( '%g', k_e ),    ...
          ',  {\itk}_{\fontsize{8}p} = ', num2str( k_p ), '  ' ];
      case this.enum.SpringModel
        imodel_info = [ this.verb.imodel{imodel},      ...
          ': {\itk}_{\fontsize{8}s} = ', ...
            sprintf( '%g', k_s ), ' N/m', ...
          ',  {\itk}_{\fontsize{8}d} = ', ...
            sprintf( '%g', k_d ), ' N m^{\fontsize{7}-1} s  ' ];
    end

    % Create caption for the symulation containing: title,
    % system info, collision model and integrator algorithm
    this.Caption = { ...
        [ '{\fontname{Arial}',                     ... % Title
          '\fontsize{10}',                         ...
          '\bf', this.Title, '}  ' ];            ...
        '{\fontsize{6} }';                         ...
        [ '\fontname{Times New Roman}',            ... % System info
          '\fontsize{11}',                         ...
          num2str( N_dim ), '-D system',           ...
          ': {\itg} = ',   num2str( this.g_n ), ' m/s^{\fontsize{7}2}', ...
          ',  {\itN}_{\fontsize{8}p} = ', num2str( N_p ) ];  ...
        '{\fontsize{2} }';                         ...
        [ '\fontname{Times New Roman}',            ... % Physical parameters 1
          '\fontsize{11}',                         ...
          '{\itC}_{\fontsize{8}D} = ',  num2str( this.C_D ), ...
          ',  {\itC}_{\fontsize{8}M} = ', num2str( this.C_M ),  ...
          ',  {\it\omega}_{\fontsize{8}M} = ', num2str( this.W_M ), ' s^{-1}', ...
          ',  {\itC}_{\it\fontsize{8}\omega} = ', num2str( this.C_W ) ]; ...
        '{\fontsize{2} }';                         ...
        [ '\fontname{Times New Roman}',            ... % Physical parameters 2
          '\fontsize{11}',                         ...
          '{\it\gamma} = ', num2str( this.gamma ),  ...
          ',  {\it\rho} = ',  num2str( this.rho ), ...
          ' kg/m^{\fontsize{7}3}', ...
          ', {\it\mu}_{\fontsize{8}part} = ',  num2str( this.mu_part ), ...
          ',  {\it\mu}_{\fontsize{8}sfc} = ', num2str( this.mu_sfc ) ];  ...
        '{\fontsize{2} }';                         ...
        imodel_info;                               ... % Interaction model
        '{\fontsize{2} }';                         ...
        [ this.verb.stepper{this.stepper},     ... % Integrator algorithm
          ': {\ith} = ',   num2str( h   ), ' s',  ...
          ',  {\itt}_{\fontsize{8}f} = ', num2str( t   ), ' s  ' ] ...
    };

end % post_UpdateCaption

%% Update animation figure and optionaly save the figure as PNG file and close
%% captured movie file

function post_CloseAnimation ()

    if ~this.Animate
        return
    end

    % Be more verbose since we don't need fast updates any more
    set( textInfo, 'String', sprintf( '{\\itN}_p = %d, {\\itt} = %.3f s', N_p, t ) );
    set( textInfo, 'Interpreter', 'tex', 'Color', 'k' );

    % Display velocity vector and particle state variables
    if this.AnnotInit && ~this.ManualNp
        for i = 1:N_p
            annotateParticle( i, struct( 'm', M(i,:), 'r', R(i,:), ...
                'x', X(i,:), 'v', V(i,:), 'w', W(i,:), 'color', FaceRgba(i,:) ) );
        end
    end

    % Remove 'cancel simulation' menu item
    delete( miCancel );
    set( mainFig, 'UserData', [] );
    drawnow;  % Flush event queue and update figure window
    set( 0, 'CurrentFigure', mainFig );  % Go back to mainFig

    % Close the movie file
    if this.CaptureMovie
        close( movieObj );            
    end

    % Annotate simulation with the caption in the top-left corner of the axes
    if ~( this.GolfTweaks && N_dim == 3 )
        axlim.x = xlim;
        axlim.y = ylim;
        caph = text( axlim.x(1) + LTick/10, axlim.y(2) - LTick/20, ...
            this.Caption, 'Interpreter', 'TeX', 'Tag', 'figCaption', ...
            'FontName', 'Helvetica', 'FontSize', 10, 'Margin', 5, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'BackgroundColor', 'none', 'EdgeColor', 'none', 'LineStyle', '-' ...
        );
        set( caph, 'UIContextMenu', get( mainFig, 'UIContextMenu' ) );
        if N_dim == 3
            axlim.z = zlim;
            set( caph, 'Position', [ axlim.x(2), LTick/5, axlim.z(2) ] );
            set( caph, 'Units', 'pixels' ); % convert to pixels i.e. fix position in 2-D
        end

        % Convert all 'Helvetica' fonts to 'Arial' and 'Times' to 'Times New Roman'
        set( findall( findall( mainFig, '-property', 'FontName' ), ...
            'FontName', 'Helvetica' ), 'FontName', 'Arial' );
        set( findall( findall( mainFig, '-property', 'FontName' ), ...
            'FontName', 'Times' ), 'FontName', 'Times New Roman' );
    end
end

%% Plot energy, linear-momentum, interactions and barycenter height

function post_PlotEnergy ()

    % Determine total number of subplots

    if track.n <= 1
        subplot_count = 0;
    else
        subplot_count = this.PlotContacts + this.PlotEnergy ...
                      + this.PlotMomentum + this.PlotHeight + this.PlotAngularW;
    end

    if subplot_count

        % Setup the figure name (windows title)
        fig_name = '';
        if this.PlotContacts
            fig_name = [ fig_name, '/Collisions' ];
        end
        if this.PlotEnergy
            fig_name = [ fig_name, '/Energies' ];
        end
        if this.PlotMomentum
            fig_name = [ fig_name, '/Linear Momentum' ];
        end
        if this.PlotHeight
            fig_name = [ fig_name, '/Height' ];
        end
        if this.PlotAngularW
            fig_name = [ fig_name, '/Angular Velocity' ];
        end
        % Get rid of leading '/' and add the simulation ID as the prefix
        fig_name = fig_name(2:end);
        fig_name = [ this.Title, ' (', this.ID, ')', ' - ', fig_name,  ]; 

        % Get the size of the screen
        screen.Rect = get( 0, 'ScreenSize' ); % gets [ left, bottom, width, height ]
        screen.L = screen.Rect(1);  screen.B = screen.Rect(2);
        screen.W = screen.Rect(3);  screen.H = screen.Rect(4);

        % Calculate position (*Warning*: Reuse screen structure of the main figure!)
        screen.fig_W = screen.H * 2/3; % Width  = 2/3 of the screen height (not width)
        screen.fig_H = screen.H * 2/3; % Height = 2/3 of the screen height
        screen.fig_L = ( screen.W - screen.fig_W ) * 0.3 + 30; % Figure left position
        screen.fig_B = ( screen.H - screen.fig_H ) * 0.6 - 30; % Figure bottom position

        % Reduce height depending of the number of the subplots
        max_plots = max( 3, subplot_count );
        plot_offset = screen.fig_H * ( max_plots - subplot_count ) / max_plots;
        screen.fig_H = screen.fig_H - plot_offset;
        screen.fig_B = screen.fig_B + plot_offset;

        % Setup paper size
        screen.pap_W = this.PlotWidth;
        screen.pap_H = this.SubplotHeight * subplot_count;

        % Create figure
        energyFig = figure( 'Name', fig_name, 'NumberTitle', 'off', ...
            'FileName', [ 'lab3_fig2_', this.ID, '.fig' ], ... % name of the FIG-file 
            'Position', [ screen.fig_L, screen.fig_B, screen.fig_W, screen.fig_H ], ...
            'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
            'PaperSize', [ screen.pap_W, screen.pap_H ], ...
            'PaperPosition', [ 0, 0, screen.pap_W, screen.pap_H ] ...
        );

        % Allocate space for individual plots axes handles
        ax = zeros( subplot_count, 1 );
        if subplot_count == 1
            ax(1) = gca;
        end
        cur_subplot = 0;  % Track the current subplot number
    end

    %% » Plot collision statistics per time-step

    if this.PlotContacts && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Collision Statistics', 'FontSize', 11 );
        set( ax( cur_subplot ) , 'YGrid', 'on' );
        ylabel( { 'Total # of contacts'; 'per time-step' }, ...
            'FontSize', 11, 'FontName', 'Times', 'Margin', 1 );

        % Calculate total number of colissions per time-step
        totColCount = sum( track.colc, 2 );

        % Plot statistics for total number of collisions
        i = find( totColCount > 0 );
        if ~isempty(i)
            stem( track.T(i), totColCount(i), 'DisplayName', 'part-part', ...
                'Color', [ 0.4 0.8 0.8 ], 'Marker', '.', ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none' );
        end

        % Plot statistics for collisions particles to surfaces
        % This plots overlays and splits previous plot into two parts: upper part with
        % particle-particle collision count and lower part with particle-surface collision
        % count.
        i = find( track.colc( :, 2 ) > 0 );
        if ~isempty(i)
            stem( track.T(i), track.colc(i,2), 'DisplayName', 'part-surface', ...
                'Color', [ 0.8 0 0 ], 'Marker', '.', ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none' );
        end

        % Finally, we have plots so we can place legend
        set( legend( 'show' ), 'Location', 'SouthEast', 'FontSize', 9 );
    end

    %% » Plot the kinetic, potential and total energy of the system over time

    if this.PlotEnergy && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Energy', 'FontSize', 11 );
        grid on;
        ylabel( '{\itE} /\fontsize{11}{\rmJ}', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Plots:
        plot( track.T, track.E(:,1), 'b-', 'LineWidth', 0.5 ); % E_k, blue, thin
        plot( track.T, track.E(:,2), 'r-', 'LineWidth', 0.5 ); % E_p, red, dot-sash, thin
        plot( track.T, track.E(:,3), 'k-', 'LineWidth', 1.5 ); % E_tot, black, thick

        % Legend:
        set( legend( '{\itE}_{\rmk}', '{\itE}_{\rmp}', '{\itE}_{\rmtot}' ), ...
            'FontSize', 10, 'FontName', 'Times', 'Location', 'SouthEast' );
    end

    %% » Plot the linear momentum of the system over time

    if this.PlotMomentum && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Linear Momentum', 'FontSize', 11 );
        grid on;
        ylabel( '{\itp} /\fontsize{11}({\rmkg \cdot m/s})', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Sum-up total linear-momentum, per time-step/dimesion (as matrix NT+1 * N_dim)
        P_tot = permute( sum( track.P, 2 ), [ 1 3 2 ] );

        % Configure plot colors and legend
        switch N_dim
            case 1, ctype = { 'z'              };
                    ltype = { 'b-'             };
            case 2, ctype = { 'x',  'z'        };
                    ltype = { 'r-', 'b-'       };
            case 3, ctype = { 'x',  'y',  'z'  };
                    ltype = { 'r-', 'g-', 'b-' };
        end

        % Plot linear-momentum components:
        for i = 1:N_dim
            plot( track.T, P_tot(:,i), ltype{i}, 'DisplayName', sprintf( ...
                '\\fontname{Times New Roman}\\fontsize{11}{\\itp}_{\\it%s\\rm,tot}', ...
                ctype{i} ...
            ));
        end
        set( legend( 'show' ), 'Location', 'SouthEast' );
    end

    %% » Plot the height of the barycenter of the system over time

    if this.PlotHeight && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end

        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Barycenter Height', 'FontSize', 11 );
        grid on;
        ylabel( '{\itz}_{\rmcm} /\fontsize{11}{\rmm}', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Find out the barycenter of the system over time as:
        %   x_cm =  Sum( x_i * m_i ) / Sum m_i

        % First extract positions and masses in given time
        X_cm = track.X;
        M_cm = repmat( permute( M, [ 3 1 2 ] ), [ track.n, 1, 1 ] );
        % Then clear masses for particles with NaN positions
        M_cm( isnan( X_cm ) ) = 0;
        X_cm( isnan( X_cm ) ) = 1;
        % Finally calculate barycenter
        X_cm = sum( M_cm .* X_cm, 2 ) ./ sum( M_cm, 2 );

        % Get the height of the barycenter and the upper envelope of the height
        height = X_cm(:,end);
        uenv = WoP.Envelope( height );

        % Plots:
        plot( track.T, height, 'b-', 'LineWidth', 0.5 ); % height: blue, thick
        plot( track.T, uenv,   'r:', 'LineWidth', 0.5 ); % envelope: red, thin
        set( legend( 'height', 'envelope' ), 'Location', 'SouthEast', 'FontSize', 9 );
        
        clearvars X_cm M_cm height uenv
    end

    %% » Plot the angular velocity 

    if this.PlotAngularW && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Angular Velocity', 'FontSize', 11 );
        grid on;
        ylabel( '{\it\omega} /\fontsize{11}({\rmrad/s})', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Configure plot colors and legend
        ctype = { 'x',  'y',  'z'  };
        ltype = { 'r-', 'g-', 'b-' };

        % Plot angular velocity components for the first particle:
        for i = 1:3
            plot( track.T, track.W(:,1,i), ltype{i}, 'DisplayName', sprintf( ...
                '\\fontname{Times New Roman}\\fontsize{11}{\\it\\omega}_{\\it%s}', ...
                ctype{i} ...
            ));
        end
        set( legend( 'show' ), 'Location', 'SouthEast' );
    end

    %% » Link axes and optionally save figure as EPS file

    if subplot_count

        % Label horizontal axis as time and link axes over time
        xlabel( '{\itt} /{\rms}', 'FontSize', 12, 'FontName', 'Times' );
        if subplot_count > 1
            linkaxes( flipdim( ax, 1 ), 'x' );
        end

        % Add context menu to export figure
        cxMenu = uicontextmenu( 'Parent', energyFig );
        uimenu( cxMenu, 'Label', 'Export as EPS', 'Callback', [ ...
            'WoP.ExportFigure( gcbf, ''-depsc2'', ', ...
            '''lab3_fig2_', this.ID, ''' );' ...
        ]);
        uimenu( cxMenu, 'Label', 'Export as PNG', 'Callback', [ ...
            'WoP.ExportFigure( gcbf, ''-dpng'', ', ...
            '''lab3_fig2_', this.ID, ''' );' ...
        ]);
        uimenu( cxMenu, 'Label', 'Export as EMF', 'Callback', [ ...
            'WoP.ExportFigure( gcbf, ''-dmeta'', ', ...
            '''lab3_fig2_', this.ID, ''' );' ...
        ]);
        set( energyFig, 'UIContextMenu', cxMenu );
        set( ax, 'UIContextMenu', cxMenu );

        % Convert all 'Helvetica' fonts to 'Arial' and 'Times' to 'Times New Roman'
        set( findall( findall( energyFig, '-property', 'FontName' ), ...
            'FontName', 'Helvetica' ), 'FontName', 'Arial' );
        set( findall( findall( energyFig, '-property', 'FontName' ), ...
            'FontName', 'Times' ), 'FontName', 'Times New Roman' );

        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig
    end

end % post_PlotEnergy


%% » Save the final values of the state variables and the traced variables
% Saves the system state variables so they could be used for initialization data 
% for another simulation.

function post_SaveData ()

    %%% Update object with the final values of the state variables
    
    this.t_f      = t;         % Simulation time
    this.t_real   = t_real;    % Elapsed real-time
    this.t_idle   = t_idle;    % Idle time
    this.t_incalc = t_incalc;  % Time spent while solving equations
    this.t_inupd  = t_inupd;   % Time spent while updating plot data
    this.t_indraw = t_indraw - t_idle - t_inupd;  % Time spent in drawing

    this.N_p = N_p;  % Number of particles (may have changed if emitting new ones)
    this.M   = M;    % Final masses
    this.R   = R;    % Final radii
    this.X   = X;    % Final positions
    this.V   = V;    % Final velocities
    this.W   = W;    % Final angular velocities
    this.FaceRgba = FaceRgba;  % Save final particle face colors

    %%% Purge time-steps that are not simulated

    if trace_T
        track.T    ( track.n+1:end )       = [];
    end
    if trace_XP
        track.X    ( track.n+1:end, :, : ) = [];
        track.P    ( track.n+1:end, :, : ) = [];
        track.W    ( track.n+1:end, :, : ) = [];
        track.colc ( track.n+1:end, : )    = [];
    end
    if trace_E
        track.E    ( track.n+1:end, : )    = [];
    end

    %%% Update traced variables

    this.track = track;   

end % function post_SaveTaces

%-------------------------------------------------------------------------------------
%% initParticle_StateVars : Initialize state variables for a new particle
% Initializes state, calculated and derived variables, and optionally a phase-space
% trajectory.
%
% Arguments:
% * i:   particle number
% * nc:  current time-step number track.n (n = 1 corresponds to t = 0)
% * pp:  structure with the particle state variables as fields:
%        'm' as mass, 'r' as radius, 'x' as position, 'v' as velocity and
%        'w' as angular velocity

function initParticle_StateVars( i, nc, pp )

    rho  = this.rho;              % Mass density of the fluid
    C_D  = this.C_D;              % Drag coefficient
    C_M  = this.C_M / this.W_M;   % Normalized Magnus coeff. See note for C_M/W_M in WoP.m
    C_W  = this.C_W;              % Angular velocity damping coefficient
    A    = pi * pp.r .^ 2;        % Cross section area of the particle
    
    %%% Initialize state, calculated and derived variables

    M   (i,:) = pp.m * ones( 1, N_dim );    % Particle mass
    R   (i,:) = pp.r;                       % Radius
    X   (i,:) = pp.x;                       % Position
    V   (i,:) = pp.v;                       % Velocity
    W   (i,:) = pp.w;                       % Angular velocity
    k_D (i,:) = 0.5 * rho * A * C_D;        % Drag coefficent
    k_M (i,:) = 0.5 * rho * A * C_M;        % Magnus force coefficient
    k_W (i,:) = 0.5 * rho * A * C_W / pp.m; % Angular velocity damping coeff.

    %%% Initialize trajectory in phase space

    if trace_XP
        track.X(nc,i,:) = X(i,:);        track.X(1:nc-1,i,:) = NaN;
        track.P(nc,i,:) = pp.m * pp.v;   track.P(1:nc-1,i,:) = 0;
        track.W(nc,i,:) = W(i,:);        track.W(1:nc-1,i,:) = 0;
    end
end

%-------------------------------------------------------------------------------------
%% initParticle_Visual : Initialize graphics object handles associated to particles
% Creates graphic handles for displaying particle and its position trajectory.
% Arguments:
% * i:  particle number

function initParticle_Visual( i )

    % Make the figure mainFig to current (i.e. target for graphics output),
    % but do not change its visibility or stacking with respect to other figures
    set( 0, 'CurrentFigure', mainFig );

    if pAsDots
        % The position of the particle in figure
        switch N_dim
            case 1,  c = [ 0, X(i) ];
            case 2,  c = X(i,:);
            case 3,  c = [ X(i,1), X(i,3) ];
        end
        % Create the graphic object, a 2D fill, representing the particle
        OH(i) = plot( c(1), c(2), ...
            '.', 'DisplayName', sprintf( 'Particle %d', i ), ...
            'Color', FaceRgba(i,1:3) .* 0.8, 'MarkerFaceColor', FaceRgba(i,1:3) ...
        );
    elseif this.PAsSpheres
        % Fix particle size
        % The position of the particle in figure
        switch N_dim
            case 1,  c = [ 0, X(i), 0 ];
            case 2,  c = [ X(i,:), 0 ];
            case 3,  c = X(i,:);
        end
        % Create the graphic object, a 3D surface, representing the particle
        % (with edge color only in 3D)
        OH(i) = surf( ...
            R(i) * US_x + c(1), R(i) * US_y + c(2), R(i) * US_z + c(3), ...
            'DisplayName', sprintf( 'Particle %d', i ), ...
            'EdgeLighting', 'gouraud', 'EdgeColor', 'none', ...
            'FaceLighting', 'gouraud', 'FaceColor', FaceRgba(i,1:3) ...
        );
        if N_dim == 3
            set( OH(i), 'EdgeColor', FaceRgba(i,1:3) );
        end
    else
        % The position of the particle in figure
        switch N_dim
            case 1,  c = [ 0, X(i) ];
            case 2,  c = X(i,:);
            case 3,  c = [ X(i,1), X(i,3) ];
        end
        % Create the graphic object, a 2D fill, representing the particle
        OH(i) = fill( R(i) * US_x + c(1), R(i) * US_y + c(2), ...
            '', 'DisplayName', sprintf( 'Particle %d', i ), ...
            'EdgeColor', FaceRgba(i,1:3) .* 0.8, ...
            'FaceColor',  FaceRgba(i,1:3), 'FaceAlpha', FaceRgba(i,4) ...
        );
    end

    % Create an empty trajectory plot object for particle (i)
    % with the same but slightly darker color as the particle
    if this.PlotPath
        trajPlot(i) = plot( NaN, NaN, '-', ...
            'DisplayName', sprintf( 'Particle %d', i ), ...
            'Color', FaceRgba(i,1:3) .* 0.8 );
        if this.GolfTweaks
            set( trajPlot(i), 'LineWidth', 1.5 );
        end
    end

end % function initParticle_Visual

%-------------------------------------------------------------------------------------
%% removeParticle : Removes particle(s) from the system
%
% Arguments:
% * ilist: list of particle numbers

function removeParticle( ilist )

    %%% Remove state and calculated variables

    M   (ilist,:) = [];    % Particle mass
    R   (ilist,:) = [];    % Radius
    X   (ilist,:) = [];    % Position
    V   (ilist,:) = [];    % Velocity
    W   (ilist,:) = [];    % Angular velocity
    k_D (ilist,:) = [];    %#ok<SETNU> Drag coefficent
    k_M (ilist,:) = [];    %#ok<SETNU> Magnus force coefficient
    k_W (ilist,:) = [];    %#ok<SETNU> Angular velocity damping coeff.

    %%% Remove trajectory in phase space

    if trace_XP
        track.X(:,ilist,:) = [];
        track.P(:,ilist,:) = [];
        track.W(:,ilist,:) = [];
    end

    %%% Remove visual representation

    if this.Animate
        removeh = OH(ilist);
        delete( removeh( ishandle( removeh ) ) );
        OH(ilist) = [];
        if this.PlotPath
            removeh = trajPlot(ilist);
            delete( removeh( ishandle( removeh ) ) );
            trajPlot(ilist) = [];
        end
    end

    %%% Update number of particles left in the system
    
    N_p = size( M, 1 );
end

%-------------------------------------------------------------------------------------
%% emitNewParticle_HiZ : Emit new particle
% Emits a new particle placing it on the top of the particle with the highest
% position.
% Arguments:
% * nc: current time-step number

function emitNewParticle_HiZ( nc )

    % Find particle with maximum height (the top-most particle)
    [~,zi] = max( X(:,end) );

    % Calculate initial state variables
    if isempty(zi)
        % If there is no the top-mostparticle, generate random position for
        % a new particle inside the box (with some margin)
        x_0 = 1.1 * this.r_p + rand( 1, N_dim ) .* ( L - 2.2 * this.r_p );
        v_0 = ZeroVec;
    else
        % Copy the state variables from the top-most particle
        x_0 = X(zi,:);
        v_0 = V(zi,:);
        % Place a new particle above the top-most particle, having the horizontal
        % initial velocity the same the top-most particle but slightly disturbed. 
        % Note that the last index 'end' always denotes the vertical coordinate.
        x_0(end) = x_0(end) + 2 * R(zi) + this.dx_emit;
        v_0(1:end-1) = v_0(1:end-1) + this.dv_emit * ( -1 + 2 * rand(1,N_dim-1) );
    end
    % Angular velocity
    switch N_dim
        case 1, w_0 = [ 0, 0, 0 ];        % no angular velicity in 1-D
        case 2, w_0 = [ 0, 0, this.w_p ]; % angular velocity around y axis
        case 3, w_0 = [ 0, this.w_p, 0 ]; % angular velocity around y axis
    end

    % Now, add particle to the system...
    N_p = N_p + 1;

    initParticle_StateVars( N_p, nc, ...
        struct( 'm', this.m_p, 'r', this.r_p, 'x', x_0, 'v', v_0, 'w', w_0 ) );

    FaceRgba(N_p,:) = this.randomFaceColor( rand ); % init new face color

    if this.Animate
        initParticle_Visual( N_p );
    end

end % function emitNewParticle

%-------------------------------------------------------------------------------------
%% addNewParticles_User : Emit new particles interactivelly
% Allow user to specify mass, radius, position and velocity of particles
% interactivelly.
% Arguments:
% * nc: current time-step number (optional; default 1)

function addNewParticles_User( nc )

    if ~nargin
        nc = 1;
    end

    % Figure header
    header = annotation( 'textbox', [ 0, 0, 1, 0.995 ], ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Margin', 5, ...
        'TextColor', [ 0.7, 0, 0 ], 'BackgroundColor', [ 0.925, 0.914, 0.847 ], ...
        'LineStyle', '-', 'EdgeColor', [ 0.5, 0.5, 0.4 ], 'FitBoxToText', 'on' );

    hue = 0;

    while true

        ip = N_p + 1;

        set( header, 'String', [ 'Click ''{\bfOK}'' to place the particle ', ...
            '{\color{gray}... or ...} Click ''{\bfCancel}'' to finish and exit' ] ...
        );

        % Ask user for particle mass and radius
        answer = inputdlg( ...
            { [ '\fontsize{10}Enter particle mass: {{\itm}_', ... % question 1
                num2str(ip), ' / kg =}' ],...
              [ '\fontsize{10}Enter particle radius: {{\itr}_', ... % question 2
                num2str(ip), ' / m =}' ], ...
              [ '\fontsize{10}Angular velocity: {{\it\omega}_', ... % question 2
                num2str(ip), ' / (rad/s) =}' ] ...
              }, ...
            [ 'Add particle ', num2str(ip) ], ... % title
            [ 1, 40; 1, 40; 1, 40 ], ... % size 2 x 40 characters
            { num2str( this.m_p ), num2str( this.r_p ), num2str( this.w_p ) }, ... % def.v
            struct( 'Resize', 'off', 'WindowsStyle', 'modal', ... % properties
                    'Interpreter', 'TeX' ) ...
        );

        if isempty( answer )  % Quit loop if user selects 'Cancel'
            break;
        end

        % Parse mass and radius
        pp.m = str2double( answer{1} ); % parse mass
        pp.r = str2double( answer{2} ); % parse radius
        pp.w = str2double( answer{3} ); % parse radius

        % Remember particle mass and radius as the default
        this.m_p = pp.m;
        this.r_p = pp.r;
        this.w_p = pp.w;

        % ... check for errors (dont' report them, just ask for a new particle)
        if isnan(pp.m) || isinf(pp.m) || pp.m <= 0
            continue
        elseif isnan(pp.r) || isinf(pp.r) || pp.r <= 0 || pp.r >= min(L)/2
            continue
        elseif isnan(pp.w) || isinf(pp.r)
            continue
        end

        % Get the particle position
        set( header, 'String', ...
            'Select the particle position {\fontname{Times New Roman}{\bfx}}' );
        drawnow; figure( mainFig ); % Update figure window pop-up our figure
        pp.x = ginput( 1 );

        % Remove the x-component of the position when in 1-D
        if N_dim == 1, pp.x(1) = []; end

        % Add a new particle to the system
        N_p = N_p + 1;
        pp.v = ZeroVec; % with zero initial velocity
        switch N_dim
            case 1, pp.w = [ 0, 0, 0 ];        % no angular velicity in 1-D
            case 2, pp.w = [ 0, 0, this.w_p ]; % angular velocity around y axis
            case 3, pp.w = [ 0, this.w_p, 0 ]; % angular velocity around y axis
        end
        pp.color = this.randomFaceColor( hue ); % and new face color
        FaceRgba(N_p,:) = pp.color;
        initParticle_StateVars( N_p, nc, pp );
        initParticle_Visual( N_p );

        % Get the the head of the velocity vector
        set( header, 'String', ...
            'Place the head of the velocity vector {\fontname{Times New Roman}{\bfv}}' );
        drawnow; figure( mainFig ); % Update figure window pop-up our figure
        pp.v = ginput( 1 );

        % Remove the x-component of the velocity when in 1-D
        if N_dim == 1, pp.v(1) = []; end

        % Set the velocity = difference between head and tail of the vector
        pp.v = pp.v - pp.x;
        
        % Reinitialize the particle state variables since we changed velocity
        initParticle_StateVars( N_p, nc, pp );

        % Display velocity vector and particle state variables
        annotateParticle( N_p, pp );

        % Get next hue from the color wheel (leaping forward among base colors)
        hue = hue + 0.31;
        if ( hue > 1 )
            hue = hue - 1;
        end
    end

    delete( header );

end % function emitNewParticle_User

%-------------------------------------------------------------------------------------
%% annotateParticle( id, pp )
% Displays velocity vector and state variables of the particle
% Arguments:
% * id: particle#
% * pp:  structure with the particle state variables as fields:
%        'm' as mass, 'r' as radius, 'x' as position, 'v' as velocity,
%        'w' as angular velocity and 'color' as face color

function annotateParticle( id, pp )

    % Calculate tail and head points for the particle velocity
    switch N_dim
        case 1
            vecX = [ 0, 0 ];
            vecY = [ pp.x(1), pp.x(1) + pp.v(1) ];
        case 2
            vecX = [ pp.x(1), pp.x(1) + pp.v(1) ];
            vecY = [ pp.x(2), pp.x(2) + pp.v(2) ];
        case 3
            return % Annotation is not supported in 3-D
    end

    % Plot the velocity vector as line 'o----.'
    color = pp.color(1:3);
    plot( vecX, vecY, ...
        '-', 'Color', color * 0.8 );
    plot( vecX(2), vecY(2), '.', ...
        'Color', color * 0.8, 'MarkerFaceColor', color );

    % Annotate vector with the the particle state variables
    idstr = num2str( id );
    text( vecX(1) - pp.r/2, vecY(1) - pp.r, {
            [ '{\itm}_', idstr, ' =', sprintf( ' %g', pp.m(1) ), ' {\rmkg}'  ]; ...
            [ '{\itr}_', idstr, ' =', sprintf( ' %g', pp.r(1) ), ' {\rmm}'   ]; ...
            [ '{\bfx}_', idstr, ' = [', sprintf( ' %.2g ', pp.x ), '] {\rmm}'   ]; ...
            [ '{\bfv}_', idstr, ' = [', sprintf( ' %.2g ', pp.v ), '] {\rmm/s}' ]; ...
            [ '{\bf\omega}_', idstr, ' = [', sprintf( ' %.2g ', pp.w ), '] {\rm1/s}' ] ...
        }, ...
        'Interpreter', 'tex', 'FontName', 'Times', 'FontSize', 9, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left' ...
    );

end % function annotateParticle

%=========================================================================================
end % method WoP.Simulate ()
