function IMAGE = transformImage_sinc(IMAGE, XIN, YIN, MATRIX)
% Transform an image by a 2x2 matrix.
% The output variable is named the same as the input variable to speed
% things up! Memory allocation is slow.

% Do a quick test to check if an identity matrix was entered.
% if so then just output the original image. If not then
% do the resampling.
if all(MATRIX ~= [1 0; 0 1])
    
    % Image dimensions
    [h, w] = size(IMAGE);
    
    % Number of points
    nPoints = numel(XIN);

    % Input points as vectors.
    xIv = reshape(XIN, 1, nPoints);
    yIv = reshape(YIN, 1, nPoints);

    % 2 x nPoints vector of coordinates.
    inputPoints = [xIv; yIv];

    % Transformed coordinates points
    interpPoints = MATRIX \ inputPoints;
    XI = interpPoints(1, :);
    YI = interpPoints(2, :);

    % Resample the image using Sinc with blackman filter. 
    % To do: update this function to accept interpolation 
    % method as an input option.
    IMAGE = reshape(sincBlackmanInterp2(double(IMAGE), XI + w/2 + 0.5, YI + h/2 + 0.5, 8, 'blackman'), h, w);
    
end % End if all(MATRIX ~= [1 0; 0 1]);


end





