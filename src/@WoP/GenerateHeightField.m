%% GenerateHeightField : Generates a random discrete terrain heigh field grid.
% Hight field is created using gausisian exponential function.
% Returns matrix in ndgrid format (x in rows, y in columns).
%
%  Usage:
%    Z = GenerateHeightField( L, tick )
%
%  Arguments:
%    ends : matrix with [ xMin, xMax; yMin, yMax ] values
%    ntck : maximum number of ticks for the shortest dimension
%    peak : peak height
%
%  Returns: 
%   Z     : a discrete random height field
%
%  Filename: @WoP/GenerateHeightField.m
%  Revision: 0.3.1
%  Date:     2012-03-27
%  Author:   Mikica B Kocic 

function Z = GenerateHeightField( ends, ntck, peak )

    % System length = first difference of ends along columns
    L = diff( ends, 1, 2 );

    % Meshgrid resolution 
    tick = min( L( L > 0 ) ) / ntck;

    % Create meshgrid for the height field
    [ X, Y ] = ndgrid( ends(1,1):tick:ends(1,2), ends(2,1):tick:ends(2,2) );

    % Create main hill at the center of field with relative size 10 ticks
    Z = gaussExp( 0.5 * L, 0 * L, 0.2 * max(L), 0, 1.0 );

    % Create small variations in the height field with normal distribution 
    % at position 0.5 +/- 0.25 and between 1 and 6 ticks in size and height 0 and 2
    % ticks.
    for i = 1:30
        roughness = -0.2 + 0.4 * rand;
        Z = Z + gaussExp( 0.5 * L, 0.25 * L, 2 * tick, 10 * tick, roughness );
    end

    % Normalize height
    zMin = min( min( Z ) );
    zMax = max( max( Z ) );
    Z = ( Z - zMin ) * peak / ( zMax - zMin );
    
    % Finally, normalize height field height
    
    %% gaussExp: generates a peak using gaussian exponential function
    %  Arguments:
    %  mu, sigma:  Peak position is normaly distributed with mean 'mu' and 
    %              standard deviation 'sigma'. 
    %  szMin, sz:  Peak size is uniformly distributed between szMin and szMin + sz.
    %  zMax:       Peak amplitude is given with zMax.

    function res = gaussExp( mu, sigma, szMin, sz, zMax ) 
        res = zMax .* exp( ...
            - ( ( X - ( mu(1) + sigma(1) * randn ) ) / ( szMin + sz * rand ) ) .^ 2 ...
            - ( ( Y - ( mu(2) + sigma(2) * randn ) ) / ( szMin + sz * rand ) ) .^ 2 ...
        );
    end

end % function GenerateHeightField
