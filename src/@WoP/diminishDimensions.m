%% Diminishes dimensions of a matrix of spatial vectors
%
% Filename: @WoP/diminishDimensions.m
% Revision: 0.1
% Date:     2012-03-19
% Author:   Mikica B Kocic

function x = diminishDimensions( x, m, n )

    % Reduce number of rows (N_p), if needed
    if m < size(x,1)
        x( m+1:end, : ) = [];
    end

    % Reduce number of columns (N_dim), if needed. 
    if n < size(x,2)
        % Note that the last column is the height, which should be kept always
        switch n
            case 1, x( :, 1:end-1 ) = []; % Remove X and opt. Y from (X,[Y,],Z)
            case 2, x( :, 2 ) = []; % Remove Y from (X,Y,Z)
        end
    end

end % function

