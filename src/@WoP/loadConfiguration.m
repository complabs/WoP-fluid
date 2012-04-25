%% Configures a simulation according to selected profile.
% Overrides default parameters depending on configuration ID.
%
%  Filename: @WoP/loadConfiguration.m
%  Revision: 0.3.3
%  Date:     2012-03-27
%  Author:   Mikica B Kocic 

function loadConfiguration( this, varargin )

    if nargin >= 2
        configId = varargin{1};
    else
        configId = 'Default';
    end

    this.Title   = [ 'Lab3: ', configId ];
    this.ID      = datestr( now, 'YYYYmmdd_HHMMSS' );
    this.imodel  = this.enum.ImpulseModel; % default for all profiles
    this.stepper = this.enum.Leapfrog;     % default for all profiles
    
    configId = lower( configId );
    switch configId

        case 'default'
            this.t_f        = 10;

        case 'demo' % 
            % Particles
            this.N_p        = 1;
            this.T_emit     = 0.25;
            this.m_p        = 1;
            this.r_p        = 0.3;
            % Model
            this.imodel     = this.enum.ImpulseModel;
            this.mu_part    = 0.0;
            this.mu_sfc     = 0.1;
            this.EnableHF   = true;
            % Stepper
            this.stepper    = this.enum.Leapfrog;
            this.t_f        = 100; % s
            this.t_inc      = 10; % s
            % Flags
            this.PAsSpheres = true;
            this.KillInactive = true;
            this.RemoveKilled = true;

        case 'golf' % 
            % System dimensions
            this.N_dim      = 3;
            this.LTick      = 10;
            this.BoxHeight  = -50;
            this.BoxWidth   = 100;
            this.BoxDepth   = 200;
            this.EnableHF   = true;
            this.GolfTweaks = true;
            % Sink
            this.R_sink     = 0.05;
            this.X_sink     = [ 50, 180, 4 ];
            % Model
            this.imodel     = this.enum.ImpulseModel;
            this.h_hard     = 0.004;
            this.IgnorePP   = true;
            this.k_e        = 0.4;
            this.k_p        = 0.2;
            this.mu_part    = 0;
            this.mu_sfc     = 0.4;
            this.gamma      = 2;
            this.C_D        = 0.3;
            this.C_M        = 0.025;
            this.C_W        = 0.3;
            this.W_M        = 63;
            % Stepper
            this.stepper    = this.enum.Leapfrog;
            this.t_f        = 40; % s
            this.t_inc      = 10; % s
            % Flags
            this.KillInactive  = true;
            this.KillOutOfBox  = true;
            this.PAsSpheres    = true;
            this.PlotPath      = true;
            this.InRealTime    = true;
            this.PctFullScreen = 0.8;
            this.MaxFigSize    = 28;
            this.PlotWidth     = 28;
            this.Silent        = true;
            % Particle
            this.m_p        = 0.045; % kg
            this.r_p        = 0.02; % m
            this.w_p        = 180; % 10800 rpm
            this.N_p        = 1;
            this.R          = this.r_p;
            this.M          = [ this.m_p, this.m_p ];
            this.X          = [ 50, 5, 5 ];
            this.V          = [ 0, 0, 0 ];
            this.W          = [ 0, 0, 0 ];
            this.FaceRgba   = [ 1, 0, 0, 1 ];

        case { 'golf-1.a', 'golf-1.b', 'golf-1.c', 'golf-1.d', 'golf-1.e' }
            % Maximum distance of the golf boll
            % System dimensions
            this.N_dim      = 2;
            this.BoxHeight  = -250;
            this.BoxWidth   = 700;
            this.LTick      = 50;
            % Model
            this.imodel     = this.enum.ImpulseModel;
            this.h_hard     = 0.01;
            this.IgnorePP   = true;
            this.k_e        = 1e-3;
            this.k_p        = 1e-2;
            this.mu_part    = 0;
            this.mu_sfc     = 1e5;
            this.EnableHF   = false;
            % Stepper
            this.stepper    = this.enum.Leapfrog;
            this.t_f        = 14; % s
            this.t_inc      = 0; % s
            % Flags
            this.KillInactive  = true;
            this.PAsSpheres    = false;
            this.PlotPath      = true;
            this.PlotAngularW  = false;
            this.InRealTime    = false;
            this.fps           = 10;
            this.PctFullScreen = 1;
            this.MaxFigSize    = 28;
            this.PlotWidth     = 28;
            % Particles
            this.m_p        = 0.045;         % Mass, kg
            this.r_p        = 0.02;          % Radius, m
            v_0             = 80;            % Initial velocity, m/s
            ang_0           = 45 * pi / 180; % Central angle,°
            ang_Range       = 10 * pi / 180; % Range ±5° around the central angle
            this.N_p        = 21;            % 0.5° between two angles
            switch configId
                case 'golf-1.a'
                    this.gamma      = 1;
                    this.C_D        = 0;
                    this.C_M        = 0;
                    this.C_W        = 0;
                case 'golf-1.b'
                    this.gamma      = 1;
                    this.C_D        = 0.3;
                    this.C_M        = 0;
                    this.C_W        = 0;
                case 'golf-1.c'
                    this.BoxHeight  = -100;
                    this.BoxWidth   = 250;
                    this.gamma      = 2;
                    this.C_D        = 0.3;
                    this.C_M        = 0;
                    this.C_W        = 0;
                    ang_0           = 37 * pi / 180; % 32-42°
                case 'golf-1.d'
                    this.BoxHeight  = -100;
                    this.BoxWidth   = 250;
                    this.gamma      = 2;
                    this.C_D        = 0.3;
                    this.C_M        = 0.025;
                    this.W_M        = 1;
                    this.C_W        = 0;
                    ang_0           = 12 * pi / 180; % 7-17°
                case 'golf-1.e'     
                    % #9 iron, 30-45°, 144 km/h, 10800 rpm
                    this.PlotAngularW = true;
                    this.LTick      = 10;
                    this.BoxHeight  = -50;
                    this.BoxWidth   = 120;
                    this.t_f        = 7; % s
                    this.gamma      = 2;
                    this.C_D        = 0.3;
                    this.C_M        = 0.025;
                    this.C_W        = 0.5;
                    this.W_M        = 40;
                    this.w_p        = 1130; % 10800 rpm
                    this.N_p        = 21; 
                    v_0             = 40; % 144 km/h (90 mph)
                    ang_0           = 25 * pi / 180;
                    ang_Range       = ang_Range * 2;
            end
            for i = 1 : this.N_p
                ang = ang_0 + ( i - 1 - this.N_p/2 ) / this.N_p * ang_Range;
                this.R(i,:)   = this.r_p;
                this.M(i,:) = [ this.m_p, this.m_p ];
                this.X(i,:) = [ this.r_p, this.r_p ];
                this.V(i,:) = [ v_0 * cos(ang), v_0 * sin(ang) ];
                this.W(i,:) = [ 0, 0, this.w_p ];
            end

        case { '1.2.a', '1.2.b' } % Energy plot, Impulse model, e = 1, Semi-implicit Euler
            % System dimensions
            this.N_dim      = 1;
            this.BoxHeight  = -5;
            this.m_p        = 1;
            this.r_p        = 0.3;
            % Particles
            this.N_p        = 1;
            this.X          = 2.3; 
            this.V          = 0;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            % Model
            this.imodel     = this.enum.ImpulseModel;
            this.k_e        = 1;
            this.k_p        = 0.01;
            this.mu_part    = 0;
            this.mu_sfc     = 0;
            this.EnableHF   = false;
            % Stepper
            this.t_f        = 50;
            this.t_inc      = 10;
            switch configId
                case '1.2.a', this.stepper = this.enum.SemiImplicitEuler;
                case '1.2.b', this.stepper = this.enum.Leapfrog;
            end
            % Flags
            this.Animate    = false;
            this.PlotEnergy = true;
            this.PlotHeight = true;

        case { '1.2.c', '1.2.d' } % Energy plot, Spring model, k_d = 0
            % System dimensions
            this.N_dim      = 1;
            this.BoxHeight  = -5;
            this.m_p        = 1;
            this.r_p        = 0.3;
            % Particles
            this.N_p        = 1;
            this.X          = 2.3; 
            this.V          = 0;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            % Model
            this.imodel = this.enum.SpringModel;
            this.k_s        = 4000;
            this.k_d        = 0;
            this.mu_part    = 0;
            this.mu_sfc     = 0;
            this.EnableHF   = false;
            % Stepper
            this.t_f        = 40;
            this.t_inc      = 10;
            switch configId
                case '1.2.c',  this.stepper = this.enum.SemiImplicitEuler;
                case '1.2.d',  this.stepper = this.enum.Leapfrog;
            end
            % Flags
            this.Animate    = false;
            this.PlotEnergy = true;
            this.PlotHeight = true;

        case { '1.3.a', '1.3.b' } % Billiard balls, Impulse model, e ? 1, Leapfrog
            % System dimensions
            this.N_dim      = 2;
            this.BoxHeight  = 2;
            this.BoxWidth   = 7;
            this.g_n        = 0;
            % Particles
            this.m_p        = 1;
            this.r_p        = 0.3;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            % 5 billiard balls in 2-D with different masses and different velocities:
            this.N_p   = 5;
            this.M = [   5,   5;  2,    2;   3,   3;   4,   4;   5,   5 ];
            this.R = [      0.3;      0.3;      0.3;      0.3;      0.3 ];
            this.X = [ 0.5, 0.3;  2,  0.3; 2.6, 0.3; 3.2, 0.3; 6.4, 0.3 ]; 
            this.V = [   1,   0;  0,    0;   0,   0;   0,   0;  -2,   0 ];
            % Model
            this.imodel     = this.enum.ImpulseModel;
            switch configId
                case '1.3.a',  this.k_e = 1.0 ; this.k_p = 0.1 ;
                case '1.3.b',  this.k_e = 0.8 ; this.k_p = 0.1 ;
            end
            this.mu_part    = 0;
            this.mu_sfc     = 0;
            this.EnableHF   = false;
            % Stepper
            this.t_f        = 5;
            this.t_inc      = 5;
            this.stepper    = this.enum.Leapfrog;
            % Flags
            this.Animate    = true;
            this.PlotContacts = true;
            this.PlotEnergy   = true;
            this.PlotMomentum = true;

        case { '2.1.a', '2.1.b' } % Billiard balls, Spring mdl, k_d?, Semi-Implicit Euler
            % System dimensions
            this.N_dim      = 2;
            this.BoxHeight  = 2;
            this.BoxWidth   = 7;
            this.g_n        = 0;
            % Particles
            this.m_p        = 1;
            this.r_p        = 0.3;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            % 5 billiard balls in 2-D with different masses and different velocities:
            this.N_p   = 5;
            this.M = [   5,   5;  2,    2;   3,   3;   4,   4;   5,   5 ];
            this.R = [      0.3;      0.3;      0.3;      0.3;      0.3 ];
            this.X = [ 0.5, 0.3;  2,  0.3; 2.6, 0.3; 3.2, 0.3; 6.4, 0.3 ]; 
            this.V = [   1,   0;  0,    0;   0,   0;   0,   0;  -2,   0 ];
            % Model
            this.imodel     = this.enum.SpringModel;
            switch configId
                case '2.1.a',  this.k_s = 4000 ; this.k_d =  0 ;
                case '2.1.b',  this.k_s = 4000 ; this.k_d = 50 ;
            end
            this.mu_part    = 0;
            this.mu_sfc     = 0;
            this.EnableHF   = false;
            % Stepper
            this.t_f        = 5;
            this.t_inc      = 5;
            this.stepper    = this.enum.SemiImplicitEuler;
            % Flags
            this.Animate      = true;
            this.PlotContacts = true;
            this.PlotEnergy   = true;
            this.PlotMomentum = true;

        case { '2.2.a', '2.2.b' } % Piling particles, Impulse model, Leapfrog, h = ?
            % System dimensions
            this.N_dim      = 1;
            this.BoxHeight  = -10;
            this.LTick      = 2;
            % Particles
            this.m_p        = 1;
            this.r_p        = 0.3;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            this.N_p        = 1;
            this.T_emit     = 1; % emit new particle every second
            this.dx_emit    = 0; % no spacing between particles
            this.M          = [ 1, 1 ];
            this.R          = 0.3;
            this.X          = 0.3;
            this.V          = 0.0;
            % Model
            this.imodel     = this.enum.ImpulseModel;
            this.k_e        = 0.9;
            this.k_p        = 0.1;
            this.mu_part    = 0;
            this.mu_sfc     = 0;
            % Stepper
            this.t_f        = 20;
            this.t_inc      = 1;
            this.stepper    = this.enum.Leapfrog;
            switch configId
                case '2.2.a',  this.h_hard = 0.01;
                case '2.2.b',  this.h_hard = 0.005;
            end
            % Flags
            this.Animate    = true;
            this.PAsSpheres = true;
            this.PlotEnergy = true;

        case { '2.3.a', '2.3.b', '2.3.c', '2.3.d', '2.3.e' }
            % Piling particles, Spring model, k_s = ?, Semi-implicit Euler
            % System dimensions
            this.N_dim      = 1;
            this.BoxHeight  = -10;
            this.LTick      = 2;
            % Particles
            this.m_p        = 1;
            this.r_p        = 0.3;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            this.N_p        = 10;
            this.T_emit     = 0.1; % emit new particle every 100 ms
            this.dx_emit    = 0; % no spacing between particles
            % Model
            this.imodel     = this.enum.SpringModel;
            switch configId
                case '2.3.a',  this.k_s =  4000 ; this.k_d = 40 ;
                case '2.3.b',  this.k_s =  8000 ; this.k_d = 40 ;
                case '2.3.c',  this.k_s =  4000 ; this.k_d =  0 ;
                case '2.3.d',  this.k_s =  8000 ; this.k_d =  0 ;
                case '2.3.e',  this.k_s = 40000 ; this.k_d = 40 ;
            end
            this.mu_part    = 0;
            this.mu_sfc     = 0;
            this.C_D        = 0;
            this.C_M        = 0;
            this.C_W        = 0;
            % Stepper
            this.t_f        = 12;
            this.t_inc      = 1;
            this.stepper    = this.enum.SemiImplicitEuler;
            % Flags
            this.Animate    = true;
            this.PAsSpheres = true;
            this.TraceVars  = true;


        case { '2.4.a', '2.4.b' } % 2D, Impulse model, stepper?
            % System dimensions
            this.N_dim      = 2;
            this.m_p        = 1;
            this.r_p        = 0.3;
            % Particles
            this.N_p        = 1;
            this.X          = [ 3, 1 ];
            this.V          = [ -0.1, 1 + rand ];
            this.T_emit     = 0.2; % emit new particle every 200 ms
            this.dv_emit    = 0.5; % variate horizontal velocity +/- 0.5 m
            % Model 
            switch configId
                case '2.4.a'
                    this.imodel = this.enum.ImpulseModel;
                    this.k_e    = 0.5;
                    this.k_p    = 0.4;
                case '2.4.b'
                    this.imodel = this.enum.SpringModel;
                    this.k_s    = 8000;
                    this.k_d    = 50;
            end
            % Stepper
            this.t_f        = 10;
            this.t_inc      = 1;
            this.stepper    = this.enum.Leapfrog;
            this.h_hard     = 0.002;
            % Flags
            this.PAsSpheres = true;
            this.MaxFigSize = 12; % Fig dimensions, in centimeters

        case { '2.5.a', '2.5.b' } % 3D, Impulse model, stepper?
            % System dimensions
            this.N_dim      = 3;
            this.m_p        = 1;
            this.r_p        = 0.3;
            % Particles
            this.N_p        = 1;
            this.X          = [ 2, 2, 1 ];
            this.V          = [ -0.1, 0.1, 1 + rand ];
            this.T_emit     = 0.2; % emit new particle every 200 ms
            this.dv_emit    = 0.5; % variate horizontal velocity +/- 0.5 m
            % Model
            switch configId
                case '2.5.a'
                    this.imodel = this.enum.ImpulseModel;
                    this.k_e    = 0.5;
                    this.k_p    = 0.4;
                case '2.5.b'
                    this.imodel = this.enum.SpringModel;
                    this.k_s    = 8000;
                    this.k_d    = 100;
            end
            % Stepper
            this.t_f        = 70;
            this.t_inc      = 1;
            this.stepper    = this.enum.Leapfrog;
            this.h_hard     = 0.002;
            % Flags
            this.MaxFigSize = 12; % Fig dimensions, in centimeters

        case 'bench'
            % System dimensions
            this.N_dim      = 3;
            this.BoxHeight  = -5;
            % Particles
            this.N_p        = 1000;
            this.X          = 2.3;
            this.V          = 0;
            % Model
            this.imodel     = this.enum.ImpulseModel;
            this.k_e        = 0.9;
            this.k_p        = 0.1;
            this.IgnorePP   = false;
            this.mu_part    = 0.01;
            this.mu_sfc     = 0.01;
            this.EnableHF   = false;
            % Stepper
            this.t_f        = 0.1;
            this.t_inc      = 0;
            this.stepper    = this.enum.Leapfrog;
            this.h_hard     = 0.001;
            % Flags
            this.Animate    = false;
            this.Silent     = true;

        otherwise
            error( 'WoP:confid', 'Unercognized configuration ID %s', configId );
    end

end % function