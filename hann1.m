function WINDOW = hann1(FILTERWIDTH, WEIGHTS, SHIFT)
% hann2 creates a 1D Hanning (cosine) Window. 
% filterDimensions = scalar specifying the width of the filter. Can be non-integer. 
%   filterDimensions = [filterHeight, filterWidth]
% weightingFactors = 1 x 2 vector specifying [weightY, weightX], where
% weightY controls the magnitude of the filter in the vertical direction,
% etc. 0 <= weightY <= 1. Weight of 0 means "don't window in this
% direction."

% % Default to even windowing.
if nargin < 3
    SHIFT = 0;
end

if nargin < 2
    WEIGHTS = [1 1];
end

% Determine dimensions and weighting factors from inputs.
filterWidth = FILTERWIDTH; % Filter height (vertical)
% filterWidth = FILTERDIMENSIONS(2); % ilter width (horizontal)
% weightY = WEIGHTINGFACTORS(1); % Vertical filter weight
% weightX = WEIGHTINGFACTORS(2); % Horizontal filter weight

shift_pix = round(SHIFT * (filterWidth - 1));

% Radial coordinate.
r = (linspace(0, pi, filterWidth/2)); % X coordinates of hann window


hannWindow1 = 1 + cos((r) .* (1/WEIGHTS(1)));
hannWindow2 = 1 + cos((r) .* (1/WEIGHTS(2)));

% hannWindow = 1 + cos(r * WEIGHTS(1));
hannWindow1(abs(r) >  pi * WEIGHTS(1))  = 0;
hannWindow2(abs(r) >  pi * WEIGHTS(2)) = 0;


% hannWindow(abs(r) > 2*pi / (2 * WEIGHTS(2))) = 0;
hannWindow = circshift(cat(2, fliplr(hannWindow1), hannWindow2), [0, shift_pix]);
r2 = (1 : length(hannWindow)) - shift_pix;
hannWindow(r2 < 0) = 0;

hannWindow(r2 > filterWidth) = 0;

% Normalize the window so its maximum value is always 1.
WINDOW = hannWindow ./ max(hannWindow(:)); % Normalize the window 
    
end