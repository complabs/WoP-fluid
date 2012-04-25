%% WoP - The World of Particles Class
% Encapsulates virtual environment framework for simulations of a system of particles.
%
% Filename: @WoP/WoP.m
% Revision: 0.1
% Date:     2012-03-19
% Author:   Mikica B Kocic

classdef WoP < handle
    
    %=====================================================================================
    %% Public properties

    properties

        %%% Simulation description

        Title   = '';
        ID      = '';
        UUID    = '';
        Caption = '';
        Description = '';

        %%% Enumueration constants, like integration method or interaction model

        enum = struct ();  % Structure holding constants
        verb = struct ();  % Structure holding verbose captions

        %%% Physical constants

        N_dim = 2;     % Number of space dimesnsions
        g_n   = 9.81;  % Standard acceleration of Gravity (intensity), m/s^2

        %%% System of particle parameters

        N_p       = 1;     % Number of particles
        N_active  = NaN;   % Number of active particles
        T_emit    = Inf;   % Emitter period (for new particles), s; may be Inf
        ManualNp  = false; % Add new particles interactivelly and quit (!)
        dx_emit   = 0.01;  % Emitter: position gap between particles, m
        dv_emit   = 0.01;  % Emitter: horizontal velocity disturbance, fraction

        m_p       = 0.045; % Particle mass, kg
        r_p       = 0.02;  % Particle radius, m
        w_p       = 10;    % Default angular velocity, rad/s

        LTick     = 1;     % Characteristic system length (distance between ticks), m
        BoxHeight = -5;    % Box width, m  (negative means: without limiting surface)
        BoxWidth  =  5;    % Box width, m
        BoxDepth  =  5;    % Box width, m
        EnableHF  = true;  % Enable hight field
        HF_z      = [];    % Height field itself (z-component)

        mu_part   = 0.1;   % Friction coefficient, in particle-particle contact
        mu_sfc    = 0.1;   % Friction coefficient, in contact with surface

        %%% Collision model physical parameters

        imodel   = NaN;   % Interaction model
        IgnorePP = false;  % Ignore particle-particle interactions

        h_soft   = 0.001; % Time-step used in linear spring model, s
        k_s      = 4000;  % Spring coefficient, N/m
        k_d      = 20;    % Spring damping factor, N/(m/s)

        h_hard   = 0.01;  % Time-step used in impulse model, s
        k_e      = 0.8;   % Impulse collision restitution (0 == max dissipation)
        k_p      = 0.3;   % Position projection fraction (1 == max projection)

        %%% Physical parameters of the fluid

        rho     = 1.3;  % The mass density of the fluid, kg/m^3
        v_fluid = [];   % Velocity of the fluid, empty vector if not used, m/s

        %%% Viscous forces parameters

        % The viscous force is f_visc = - k_D * |V|^gamma * V/|V|
        % where the viscous damping coefficient k_D = 1/2 * C_D * rho * A
        % and A is the cross-section area of an object (A = pi * r^2 for a sphere).
        % C_D depends on geometrical shape and surface material.

        gamma = 1.0;  % Type of the flow: 1 == laminar, 2 == turbulent
        C_D   = 0.3;  % Drag coefficient
        
        % The rotational drag angular acceleration causing angular velocity decay is
        % w_dot = - k_W * ( |V| + R *|W| )^(gamma-1) * W, where the damping coefficient
        % k_W = 0.5 * C_W * rho * A / m, and A is the cross-section are of an 
        % object, C_W is non-dimensional constants that depends on geometrical shape 
        % (incl. object's moment of inertia) and surface material, and k_W 
        % angular damping coefficient is given in m^-1.

        C_W  = 0.1;  % Angular velocity damping coefficient

        %%% Magnus force parameters

        % The Magnus force is f_M = k_M * |V| * ( W × V )
        % where k_M = 1/2 * rho * A * C_M / |W_M|, rho is mass density of the fulid,
        % A is the cross-section area, C_M is coefficient depending on the geomatry
        % (and other non included parameters) and |W_ref| is the referent angular 
        % velocity intensity for given C_M. This means that Magnus force is has a 
        % almost linear dependence of the angular velocity around the central intensity 
        % |W_M| for which C_M is given.

        C_M  = 0.025;  % The magnus force coefficient, non-dimensional value
        W_M  = 10;     % Referent angular velocity intensity for given C_M, rad/s

        %%% Results control flags

        TraceVars     = false ; % Trace trajectory and energy levels data
        PlotContacts  = false ; % Plot number of collisions (contacts) per time-step
        PlotEnergy    = false ; % Plot kinetic, potential and total energy over time
        PlotHeight    = false ; % Plot envelope of the first particle
        PlotMomentum  = false ; % Plot linear momentum of the system
        PlotAngularW  = false ; % Plot angular velocity of the first particle

        %%% Visualization control flags

        Animate      = true  ; % Animate particles
        InRealTime   = true  ; % Try animation in real-time, if possible
        PlotPath     = false ; % Plot particle trajectories
        PAsSpheres   = false ; % Display particles as spheres; othewise circles
        Wiggle3D     = false ; % Use 'Wiggle Stereoscopy' to convay 3-D depth
        AnnotInit    = false ; % Annotate the initial particle state variables
        AnnotFinal   = false ; % Annotate the final particle state variables
        Silent       = false ; % Do not print debugging info during simulation
        KillInactive = false ; % Make non-moving particles near the ground inactive
        KillOutOfBox = false ; % Make particles inactive if they go out of bounds
        RemoveKilled = false ; % Remove inactive particles 
        GolfTweaks   = false ; % Golf adapted setup, e.g. thicker lines etc.
        
        %%% Particle sink

        X_sink  = [];  % Sink position (only x, y coordinates; z is ignored)
        R_sink  = [];  % Sink radius

        %%% Animation parameters

        % Valid combinations of mainFig/mainAxes are:
        % 1) mainAxes and textInfo in existing figure or gui panel
        % 2) mainFig handle of existing figure
        % A new mainFig is created if mainFig/mainAxes are not specified
        % and the new mainFig is reused across simulations.

        mainFig;     % Main figure handle
        mainAxes;    % Main axes handle
        textInfo;    % Handle of a text object used to display animation progress

        fps   = 25;  % Animation frame rate, Hz
        N_usf = 24;  % Number of faces on a sphere
        N_uck = 48;  % Number of segments in a circle

        PctFullScreen = 0.5;  % Relation of the animation figure size to screen size
        
        %%% Saving figures and capturing animation flags

        CaptureMovie  = false ; % Capture animation as movie
        MaxFigSize    = 17    ; % Max width or height of saved figure, in centimeters
        PlotWidth     = 17    ; % Width of the plot, in centimeters
        SubplotHeight = 5     ; % Height of the subplot (when stacked), in centimeters

        %%% The state variables (by default empty matrices)

        M;          % Masses:       [ N_p * N_dim double ]
        R;          % Radii:        [ N_p * 1     double ]  NaN for dead particles 
        X;          % Positions:    [ N_p * N_dim double ]
        V;          % Velocities:   [ N_p * N_dim double ]
        W;          % Angular vel.: [ N_p * 3     double ]  always 3D (!)
        FaceRgba;   % Face color:   [ N_p * 4     double ]  RGB + alpha

        %%% Integration parameters

        Running  = false; % True while the simulation is running
        stepper  = NaN;   % Integrator method
        h        = Inf;   % Integrator time-step, s
        t_0      = 0;     % Initial time, s
        t_f      = Inf;   % Final time, s
        t_inc    = 0;     % Final time increment, when continuing simulation, s
        t_real   = NaN;   % Elapsed real time, s
        t_idle   = NaN;   % Paused (when idle) time, s
        t_incalc = NaN;   % Time spent while solving equations, s
        t_indraw = NaN;   % Time spent while rendering animation, s
        t_inupd  = NaN;   % Time spent wwhile updating plot data, s

        %%% The structure with the trace of the system variables (over all simulations)

        % Fields:
        % track.n     : Number of samples
        % track.T     : Elapsed time
        % track.X     : Position trajectory in phase-space
        % track.P     : Linear momentum trajectory in phase-space
        % track.colc  : Number of collisions per time-step where:
        %               1st column contains particle-to-particle collision count and 
        %               2nd column contains particle-to-surface collision count
        % track.E_k   : Kinetic energy
        % track.E_p   : Potential energy
        % track.E_tot : Total energy
        
        track = struct( 'n', 0 );   % Empty at the beginning

    end % public properties

    %=====================================================================================
    %% Read-only public properties
    
    properties ( SetAccess = protected, GetAccess = public )
        
        %%% The initial state variables (saved at t_0 by Simulate())
        M_0;         % Masses:       [ N_p * N_dim double ]
        R_0;         % Radii:        [ N_p * 1     double ]
        X_0;         % Positions:    [ N_p * N_dim double ]
        V_0;         % Velocities:   [ N_p * N_dim double ]
        W_0;         % Angular vel.: [ N_p * N_dim double ]
        FaceRgba_0;  % Face color:   [ N_p * 4     double ]

    end % Read-only public properties

    %=====================================================================================
    %% Constructor

    methods
        function this = WoP( varargin )

            if usejava('jvm')
                this.UUID = char( java.util.UUID.randomUUID () );
            end

            %%% Initialize enumeration constants in structure (approach with structucre is 
            %%% much faster than having classdef with enumerations)

            % Interaction model enumeration
            this.enum.SpringModel  = 1;   % Linear-spring interaction model
            this.enum.ImpulseModel = 2;   % Impulse collision model

            % Integrator method enumeration
            this.enum.ForwardEuler      = 1;   % Forward Euler integrator
            this.enum.SemiImplicitEuler = 2;   % Semi-implicit Euler integrator
            this.enum.Leapfrog          = 3;   % Leapfrog integrator

            % Verbose captions for enumeration constants
            this.verb.imodel  = { 'Linear-spring model'; 'Impulse collision model' };
            this.verb.stepper = { 'Forward Euler'; 'Semi-implicit Euler'; 'Leapfrog' };

            %%% Customize configuration depending on the selected profile

            this.loadConfiguration( varargin{:} );

            this.addlistener( 'ObjectBeingDestroyed', ...
                @(src,evtdata)disp([ 'ObjectBeingDestroyed: ', src.ID ]) ...
            );            
        end
    end

    %=====================================================================================
    %% Inline public methods

    methods ( Access = public )

        % delete(): The class destructor
        function delete( this )
            StopSimulation( this );
            disp([ 'Purged ', this.Title, ', ID: ', this.ID ]);
        end

        % StopSimulation(): Breaks current main simulation loop
        function StopSimulation( this )
            this.Running = false; % The simulation will run while this flag is true
        end
    end

    %=====================================================================================
    %% Public methods
    
    methods ( Access = public )
        Simulate( this )
        ResetSimulation( this )
        % varargout = FastSim( this, varargin )
    end

    %=====================================================================================
    %% Public static methods

    methods ( Static, Access = public )
        configIds = GetConfigIDs ()
        varargout = Envelope( sig, threshold )
        ExportAnnotation( figh, tag, dtype, filename ) 
        ExportFigure( figh, dtype, filename )
        [ X, Y, HF ] = GenerateHeightField( L, tick, peak )
    end

    %=====================================================================================
    %% Private methods

    methods ( Access = private )
        color = randomFaceColor( this, hue )
        loadConfiguration( this, varargin )
    end

    %=====================================================================================
    %% Private static methods

    methods ( Static, Access = private )
        fixPsFonts( filename )
        spatialVec = diminishDimensions( spatialVec, m, n )
    end

end % classdef WoP

