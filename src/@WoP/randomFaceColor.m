%% Generates a random color to be used as particle face color
% Returns a vector with 4 components: the first 3 components are the RGB color components 
% and the 4-th component is alpha (transparency)
%
% Filename: @WoP/randomFaceColor.m
% Revision: 0.2
% Date:     2012-03-19
% Author:   Mikica B Kocic 

function color = randomFaceColor( this, hue )

    if this.PAsSpheres
        % random hue with max value, but reduced saturation
        color = [ hue, 0.4 + 0.3 * rand, 1 ];
        alpha = 1;
    else
        % random hue with reduced value and saturation (pastell colors)
        color = [ hue, 0.5 + 0.3 * rand, 0.7 + 0.3 * rand ]; 
        alpha = 0.8;
    end

    % Return RGB + alpha
    color =  [ hsv2rgb( color ), alpha ];

end % function