function IMAGEOUT = deformImage(IMAGEIN, X, Y, U, V)

% Image dimensions
[imageHeight, imageWidth] = size(IMAGEIN); 

% Create the pixel coordinates.
[xi_integer, yi_integer] = meshgrid(1:imageWidth, 1:imageHeight);

% Shift the pixel coordinates by 0.5 pixels
XI = xi_integer - 0.5;
YI = yi_integer - 0.5;

% Create interpolation structures for the velocity field.
interpolant_tx = griddedInterpolant(Y, X, U, 'spline', 'linear');
interpolant_ty = griddedInterpolant(Y, X, V, 'spline', 'linear');

% This is the velocity field upsampled to every pixel.
UI = interpolant_tx(YI, XI);
VI = interpolant_ty(YI, XI);

% These are the coordinates at which to resample image 2.
XD = XI + UI;
YD = YI + VI;

% Resample the images
IMAGEOUT = sincBlackmanInterp2(IMAGEIN, XD + 0.5, YD + 0.5, 8, 'blackman');

end