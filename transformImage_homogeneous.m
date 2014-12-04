function IMAGE = transformImage_homogeneous(IMAGE, XIN, YIN, MATRIX, COMPILED)
% Transform an image by a 3x3 homogeneous matrix.
% The output variable is named the same as the input variable to speed
% things up! Memory allocation is slow.

% Do a quick test to check if an identity matrix was entered.
% if so then just output the original image. If not then
% do the resampling.
if ~all(MATRIX(:) == [1 0 0 0 1 0 0 0 1]')
    
    % Image dimensions
    [h, w] = size(IMAGE);
    
    % Number of points
    nPoints = numel(XIN);

    % Input points as vectors.
    xIv = reshape(XIN, 1, nPoints);
    yIv = reshape(YIN, 1, nPoints);

    % 3 x nPoints vector of coordinates.
    inputPoints = [xIv; yIv; ones(1, nPoints)];

    % Transformed coordinates points
    interpPoints = MATRIX \ inputPoints;
    XI = interpPoints(1, :);
    YI = interpPoints(2, :);

    % Resample the image using either interp2 (matlab) or the compiled
    % interpolator
    % ba_interp2
    if COMPILED
        IMAGE = (reshape(ba_interp2(double(IMAGE), XI + w/2 + 0.5, YI + h/2 + 0.5, 'cubic'), h, w))';
    else
        IMAGE = reshape(interp2(XIN, YIN, double(IMAGE), XI, YI, 'bicubic', 0), h, w);
    end
    
end % End if all(MATRIX ~= [1 0; 0 1]);


end





