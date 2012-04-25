%% Finds upper and lower envelopes of a given signal
%
% Filename: @WoP/Envelope.m
% Revision: 0.3
% Date:     2012-03-19
% Author:   Mikica B Kocic
%
% Usage:
%   upperenv = envelope( sig )
%   upperenv = envelope( sig, threshold )
%   [upperenv,lowerenv] = envelope( sig )
%   [upperenv,lowerenv] = envelope( sig, threshold )
%
% Arguments:
% * sig    : vector of input signal (either as 1xN or Nx1 matrix)
% * threshold : ignore signal variations bellow this value; default 1e-4
%
% Returns:
% * upperenv : upper envelope of the input signal
% * lowerenv : lower envelope of the input signal

function varargout = Envelope( sig, threshold )

    if nargin < 2
        threshold = 1e-4;
    end

    % Find out the first derivative of the signal
    delta = diff( sig );

    % Flatten signal variations bellow the threshold
    delta( abs(delta) < threshold ) = 0;
    delta = sign( delta );

    % Remove all regions where signal remains constant
    for xi = 2 : length( delta )
        if delta( xi ) == 0
            delta( xi ) = delta( xi - 1 );
        end
    end

    % Find out the second derivative
    delta = diff( delta );

    % Determine local maximum and minimum points
    upper_ind = find( delta < 0 ) + 1;  % maximum if f''(x) < 0
    lower_ind = find( delta > 0 ) + 1;  % minimum if f''(x) > 0

    first = 1;
    last = length(sig);
    xi = first : last;

    if length(xi) <= 1
        if nargout >= 1
            varargout{1} = sig;
        end
        if nargout >= 2
            varargout{2} = sig;
        end
        return
    end

    if length( upper_ind ) < 2
        upper_ind = [ first, upper_ind(:), last ];
    end
    if length( lower_ind ) < 2
        lower_ind = [ first, lower_ind(:), last ];
    end

    if nargout >= 1
        varargout{1} = interp1( upper_ind, sig(upper_ind), xi, 'linear', 'extrap' );
    end
    if nargout >= 2
        varargout{2} = interp1( lower_ind, sig(lower_ind), xi, 'linear', 'extrap' );
    end

end % function