function [x_grid_rect, y_grid_rect, x_grid_01, y_grid_01, x_grid_02, y_grid_02] = discreteWindowOffset(X_OLD, Y_OLD, U, V, ProcessingParameters) 
% Inputs: 
%   XGRID and YGRID are the matrices or vectors of the column and row 
%   grid points prior to shifting. These Should be in the format of meshgrid.
% 
%   U and V are the velocity fields by which to shift the grid points.
%   The size of U and V need not be the same as the size of XGRID and
%   YGRID. If the sizes of U and V are equal to the sizes of XGRID and
%   YGRID, then the grid will be shifted by exactly U and V; if those
%   sizes are not equal, then the shifts applied to XGRID and YGRID 
%   will be interpolated (linear) from U and V.  
%
% OUTPUTS
%   X_GRID_01 and Y_GRID_01 are the shifted grids for the first image in
%   the pair.
%   
%   X_GRID_02 and Y_GRID_02 are the shifted grids for the second image in
%   the pair.
%   

% Specify the DWO differencing scheme (forward, central, or backward)
DwoMethod = ProcessingParameters.FmcDifferenceMethod;

% Read image dimensions
imageHeight = ProcessingParameters.Images.Height;
imageWidth  = ProcessingParameters.Images.Width;

% Interrogation Region height and width
regionHeight = ProcessingParameters.InterrogationRegion.Height;
regionWidth  = ProcessingParameters.InterrogationRegion.Width;

% Grid spacing parameters
gridSpacingY = ProcessingParameters.Grid.Spacing.Y;
gridSpacingX = ProcessingParameters.Grid.Spacing.X;

% Make sure the grid buffer is at least half the size of the interrogation region
gridBufferY = max(ProcessingParameters.Grid.Buffer.Y, ceil(regionHeight / 2));
gridBufferX = max(ProcessingParameters.Grid.Buffer.X, ceil(regionWidth / 2));

% Create the initial rectangular grid
[x_grid_rect, y_grid_rect] = gridImage([imageHeight, imageWidth], [gridSpacingY gridSpacingX], gridBufferY, gridBufferX);

% Interpolate the shift-field onto the input grid. U and V are the
% non-integer horizontal and vertical displacements by whose rounded values
% the grid is shifted.
% "Spline" is the interpolation method and
% "Linear" is the extrapolation method.
% We are using griddedInterpolant because interp2 doesn't have a
% good extrapolation method input (just scalar values). 
interpolant_U = griddedInterpolant(Y_OLD,X_OLD, U, 'spline','linear');
interpolant_V = griddedInterpolant(Y_OLD,X_OLD, V, 'spline','linear');
gridShiftX = interpolant_U(y_grid_rect, x_grid_rect);
gridShiftY = interpolant_V(y_grid_rect, x_grid_rect);

% Create logical expressions to specify the differencing scheme
isBackwardDifference = ~isempty(regexpi(DwoMethod, 'ba'));
isForwardDifference  = ~isempty(regexpi(DwoMethod, 'fo'));

% Shift the grid coordinates
if isBackwardDifference
    % In case of backward differencing...
    
    % Shift the grid points from the first image by the negative
    % of the displacement field
    x_grid_01 = x_grid_rect -1 * round(gridShiftX);
    y_grid_01 = y_grid_rect -1 * round(gridShiftY);

    % Keep the original grid points from the second image.
    x_grid_02 = x_grid_rect;
    y_grid_02 = y_grid_rect;
    
elseif isForwardDifference
    % In case of forward differencing...
    
    % Keep the original grid points from the first image.
    x_grid_01 = x_grid_rect;
    y_grid_01 = y_grid_rect;
    
    % Shift the grid points from the second image by the input displacement
    % field.
    x_grid_02 = x_grid_rect + round(gridShiftX);
    y_grid_02 = y_grid_rect + round(gridShiftY);
    
else
    % In (default) case of central differencing... 
    
    % Shift the grid points from the first image by -1/2 times the input
    % displacement field.
    x_grid_01 = x_grid_rect - round(1/2 * gridShiftX);
    y_grid_01 = y_grid_rect - round(1/2 * gridShiftY);
    
    % Shift the grid points from the second image by +1/2 times the input
    % displacement field
    x_grid_02 = x_grid_rect + round(1/2 * gridShiftX);
    y_grid_02 = y_grid_rect + round(1/2 * gridShiftY);

end

% These lines prevent the grid shifts from placing grid points that
% would result in regions that extend past the image.
x_grid_01((x_grid_01 - floor(regionWidth/2)) < 1 ) = floor(regionWidth/2);
x_grid_01((x_grid_01 +  ceil(regionWidth/2)) > imageWidth) = imageWidth - ceil(regionWidth/2);

x_grid_02((x_grid_02 - floor(regionWidth/2)) < 1 ) = floor(regionWidth/2);
x_grid_02((x_grid_02 +  ceil(regionWidth/2)) > imageWidth) = imageWidth - ceil(regionWidth/2);

y_grid_01((y_grid_01 - floor(regionHeight/2)) < 1) = floor(regionHeight/2);
y_grid_01((y_grid_01 +  ceil(regionHeight/2)) > imageHeight) = imageHeight - ceil(regionHeight/2);

y_grid_02((y_grid_02 - floor(regionHeight/2)) < 1) = floor(regionHeight/2);
y_grid_02((y_grid_02 +  ceil(regionHeight/2)) > imageHeight) = imageHeight - ceil(regionHeight/2);

end