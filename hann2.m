function hannWindow = hann2(FILTERDIMENSIONS, WEIGHT)
% hann2 creates a 2D Hanning (cosine) Window. 
% filterDimensions = 1 x 2 vector specifying the height and width of the filter. Can be non-integer. 
%   filterDimensions = [filterHeight, filterWidth]
% weightingFactors = 1 x 2 vector specifying [weightY, weightX], where
% weightY controls the magnitude of the filter in the vertical direction,
% etc. 0 <= weightY <= 1. Weight of 0 means "don't window in this
% direction."

% % Default to even windowing.
% if nargin < 2
%     WEIGHT = [1 1];
% end

% Determine dimensions and weighting factors from inputs.
filterHeight = FILTERDIMENSIONS(1); % Filter height (vertical)
filterWidth = FILTERDIMENSIONS(2); % ilter width (horizontal)
% weightY = WEIGHTINGFACTORS(1); % Vertical filter weight
% weightX = WEIGHTINGFACTORS(2); % Horizontal filter weight

% Create vectors that span 0 : 3 pi of length filterHeight and filterWidth, at which the weighted cosine is
% evaluated.
hannCoordinatesY = (linspace(-pi, pi, filterHeight))'; % Y coordinates of hann window
hannCoordinatesX = (linspace(-pi, pi, filterWidth)); % X coordinates of hann window

% Create a mesh grid of coordinates.
[hannCoordinatesX, hannCoordinatesY] = meshgrid(hannCoordinatesX, hannCoordinatesY); % Make 2-D coordinates

[~, r] = cart2pol(hannCoordinatesX, hannCoordinatesY);


hannWindow = 1 + cos(r * WEIGHT);
hannWindow(r > 2*pi / (2 * WEIGHT)) = 0;


% % Evaluate the cosine function at the Hanning coordinates.
% % hannWindow =  (1 + weightY * cos(hannCoordinatesY)) * (1 + weightX * cos(hannCoordinatesX)); % Calculate window
% hannWindow =  1 + (weightY * cos(hannCoordinatesY)) * (weightX * cos(hannCoordinatesX));

% Normalize the window so its maximum value is always 1.
hannWindow = hannWindow ./ max(hannWindow(:)); % Normalize the window 
    
end