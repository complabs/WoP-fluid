%% Resets initial time back to t_0 = 0 and removes final state variables.
%
% Filename: @WoP/ResetSimulation.m
% Revision: 0.3
% Date:     2012-03-20
% Author:   Mikica B Kocic 

function ResetSimulation( this )

    this.t_0 = 0; % Reset initial time

    % Clear the track of the traced variables
    
    this.track = struct( 'n', 0 );

    % Reload the final state variables from the remembered initial conditions

    % Masses
    this.M = this.M_0;
    N_p = min( this.N_p, size( this.M_0, 1 ) );

    % Radii
    this.R = this.R_0;
    N_p = min( N_p, size( this.R_0, 1 ) );

    % Positions
    this.X = this.X_0;
    N_p = min( N_p, size( this.X_0, 1 ) );

    % Velocities
    this.V = this.V_0;
    N_p = min( N_p, size( this.V_0, 1 ) );

    % Angular velocities
    this.W = this.W_0;
    N_p = min( N_p, size( this.W_0, 1 ) );

    % Face colors
    this.FaceRgba = this.FaceRgba_0;

    % Set this.N_p to the minimal number of particles found in initial state vars
    if N_p > 0
        this.N_p = N_p;
    end
    
end % function