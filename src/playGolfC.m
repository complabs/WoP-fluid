%% The Golf Game - Command Window Application
%
%  Filename: @WoP/playGolf.m
%  Revision: 0.2
%  Date:     2012-03-31
%  Author:   Mikica B Kocic

function varargout = playGolfC

    close all;            % Close all figures

    sim = WoP( 'golf' );  % Create new simulation

    sim.t_f = 20;         % Setup max final time
    sim.Simulate;         % Place initial ball on the ground

    while true

        sim.R = sim.r_p; % Reset NaN radius

        % Add random wind

        sim.v_fluid  = [ 4 + randn * 2, -2 + randn * 0.5, 0 ]; 

        % Display wind and position status

        disp( '----------------------------------------------------------' );
        fprintf( 'Wind = ' ); fprintf( ' %9g', sim.v_fluid ); fprintf( '\n' );
        fprintf( 'X    = ' ); fprintf( ' %9g', sim.X       ); fprintf( '\n' );
        disp( '----------------------------------------------------------' );
        disp( 'Note: 0 degrees is along depth (Y-axis direction)' )
        disp( '----------------------------------------------------------' );

        % Tak inputs from user

        dir  = input( 'Direction (degrees) ? ' );
        v0   = input( 'Velocity  (km/h)    ? ' );
        loft = input( 'Loft      (degrees) ? ' );
        rpm  = input( 'Spin      (rpm)     ? ' );

        % Convert units

        v0    = v0   / 3.6;     %   1 m/s = 3.6 km/h
        w0    = rpm  * 2*pi/60; %  60 rpm = 2*pi rad/s 
        loft  = loft * pi/180;  % 180 deg = pi rad
        dir   = dir  * pi/180;  % 180 deg = pi rad

        % Setup initial state variables

        sim.V = [ v0 * sin(dir) * cos(loft) ,  ... % Vx
                  v0 * cos(dir) * cos(loft) ,  ... % Vy
                  v0 * sin(loft)            ]; ... % Vz

        sim.W = [ w0, 0, 0 ];

        sim.t_f = sim.t_f + 20;  % Advance max. final time

        disp( '----------------------------------------------------------' );
        fprintf( 'V    = ' ); fprintf( ' %9g', sim.V ); fprintf( '\n' );
        fprintf( 'W    = ' ); fprintf( ' %9g', sim.W ); fprintf( '\n' );
        disp( '----------------------------------------------------------' );

        sim.Simulate;

        % In case of collision with box surfaces, quit
        if sim.track.colc(end,2)
            fprintf( '*** Out of Bounds ***\n\n' )
            break
        elseif sim.track.colc(end,3)
            fprintf( '*** IN HOLE *** Congratulations!\n\n' )
            break
        end

    end

    if nargout >= 1
        varargout{1} = sim;
    end

end
