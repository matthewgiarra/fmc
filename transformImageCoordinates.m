function [YOUT, XOUT] = transformImageCoordinates(TRANSFORM, XGRID, YGRID, METHOD, CENTER)
% TRANSFORM = affine matrix
% IMAGESIZE = size of image (pixels); IMAGESIZE = [nRows nCols]

% Default to not cropping the image.
if nargin < 4
    METHOD = 'full';
end

% Count rows and columns in the original image 
[nr, nc] = size(XGRID);

% Determine the center of the transformation. Default to center of grid
if nargin < 5
    xc = (nc + 1) / 2;
    yc = (nr + 1) / 2;
else
    yc = CENTER(1);
    xc = CENTER(2);
end

% Shift grid origin to the center of grid
xPoints = XGRID - xc;
yPoints = YGRID - yc;

% Reshape coordinate matrices into vectors
xPointsVect = reshape(xPoints, 1, numel(xPoints));
yPointsVect = reshape(yPoints, 1, numel(yPoints));

% Apply the transformation matrix to the shifted points
transformedPoints = TRANSFORM * [xPointsVect; yPointsVect; ones( 1, length(xPointsVect ) ) ];

% Extract the X- and Y coordinates of the transformed points.
Xt = transformedPoints(1, :);
Yt = transformedPoints(2, :);

% Turn transformed X- and Y- coordinates back into matrices and shift the origin back to (1, 1)
x = reshape(Xt, nr, nc) + xc;
y = reshape(Yt, nr, nc) + yc;

% Crop the image if specified.
if strcmp(METHOD, 'crop')
    XOUT = (x(x >= min(XGRID(:)) & x <= max(XGRID(:)) & y >= min(YGRID(:)) & y <= max(YGRID(:))));
    YOUT = (y(x >= min(XGRID(:)) & x <= max(XGRID(:)) & y >= min(YGRID(:)) & y <= max(YGRID(:))));
    
else % Don't crop the image if cropping wasn't specified.
    XOUT = x;
    YOUT = y;
    
end


end



