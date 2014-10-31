
function [REGION_MATRIX_01, REGION_MATRIX_02, GRID_SHIFT_X, GRID_SHIFT_Y] = discreteWindowOffset_prana(X_02, Y_02, X_01, Y_01, U_01, V_01, IMAGE_01, IMAGE_02, REGIONWIDTH, REGIONHEIGHT)
% [regionMatrix_01, regionMatrix_02] = discreteWindowOffset(X, Y, U, V, Image1, Image2, regionWidth, regionHeight)
% This function performs central-difference discrete window offset for PIV images.
% The inputs are the original grid points, velocity estimates, PIV images, and
% dimensions of the interrogation regions. The outputs are two M x N x P
% matrices where M and N are the height and width of the new interrogation
% regions and P is the number of new interrogation regions.
%
% Grid points in the first (of two) image are shifted by negative one half
% of the input velocity field, and grid points in the second image are
% shifted by positive one half of the input velocity field. 
%
% INPUTS
%   X = Matrix specifying the integer column coordinates of the grid
%   points prior to shifting.
%
%   Y = Matrix specifying the integer row coordinates of the grid
%   points prior to shifting.
%
%   U = Vector of matrix specifying the horizontal components of ths
%   velocity field on which the discrete window offset is based.
%
%   V = Vector of matrix specifying the vertical components of ths
%   velocity field on which the discrete window offset is based.
%   
%   IMAGE_01 = Matrix containing the image data for the first image.
%
%   IMAGE_02 = Matrix containing the image data for the second image.
%
%   REGIONWIDTH = Integer width of the new interrogation regions.
%
%   REGIONHEIGHT = Integer height of the new interrogation regions.
%
% OUTPUTS
%   REGION_MATRIX_01 and REGION_MATRIX_02 are matrices of all of the new
%   interrogation regions from the first and second images, respectively. 
%   The dimensions of these matrices are [REGIONHEIGHT, REGIONWIDTH, N]
%   where N is the number of interrogation regions (i.e. grid points). 
%
% Variable names of inputs and outputs of this function are capitalized.

% Determine the image dimensions
[imageHeight, imageWidth] = size(IMAGE_01);

% Interpolate the original velocity field onto the new grid to determine
% the pixel shifts to apply to each grid point.
GRID_SHIFT_X = interp2(X_01, Y_01, U_01, X_02, Y_02, 'spline', 0);
GRID_SHIFT_Y = 1 * interp2(X_01, Y_01, V_01, X_02, Y_02, 'spline', 0);

% Convert the grid shift matrices to vectors for easy looping.
gridShiftX = GRID_SHIFT_X(:);
gridShiftY = GRID_SHIFT_Y(:);

% Change the input coordinate and velocity matrices to vectors for easy
% loopinbg. 
x = X_02(:);
y = Y_02(:);

% Determine the number of grid points.
nGridPoints = length(x);

% Initialize matrices to hold the shifted interrogation regions.
REGION_MATRIX_01 = zeros(REGIONHEIGHT, REGIONWIDTH, nGridPoints);
REGION_MATRIX_02 = zeros(REGIONHEIGHT, REGIONWIDTH, nGridPoints);

% Calculate the integer coordinates of the shifted grid points
% in the first image. 
x_shift_01 = x - floor(round(gridShiftX)/2);
y_shift_01 = y - floor(round(gridShiftY)/2);

% Calculate the integer coordinates of the shifted grid points
% in the second image. 
x_shift_02 = x + ceil(round(gridShiftY)/2);
y_shift_02 = y + ceil(round(gridShiftY)/2);

% Left and ridge edges of the new regions in the first image.
xmin1 = x_shift_01 - ceil(REGIONWIDTH / 2) + 1;
xmax1 = x_shift_01 + floor(REGIONWIDTH / 2);

% Top and bottom edges of the new regions in the first image.
ymin1 = y_shift_01 - ceil(REGIONHEIGHT / 2) + 1;
ymax1 = y_shift_01 + floor(REGIONHEIGHT / 2);

% Left and right edges of the new regions in the second image.
xmin2 = x_shift_02 - ceil(REGIONWIDTH / 2) + 1;
xmax2 = x_shift_02 + floor(REGIONWIDTH / 2);

% Top and bottom edges of the new regions in the second image.
ymin2 = y_shift_02 - ceil(REGIONHEIGHT / 2) + 1;
ymax2 = y_shift_02 + floor(REGIONHEIGHT / 2);

% This extracts the image intensity values at the locations corresponding
% to the new windows. 
for n = 1 : nGridPoints

    % Extract the interrogation regions from the images.
    shiftedSubRegion_01 = IMAGE_01( max(1, ymin1(n)) : min(imageHeight, ymax1(n)), max(1, xmin1(n)) : min(imageWidth, xmax1(n)));
    shiftedSubRegion_02 = IMAGE_02( max(1, ymin2(n)) : min(imageHeight, ymax2(n)), max(1, xmin2(n)) : min(imageWidth, xmax2(n)));
    
    % If the dimensions of the first shifted region are less than the specified
    % dimensions of the interrogation regions, then Pad the image subregion
    % with zeros so that its final dimensions are consistent with those
    % specified.
    if all([REGIONHEIGHT, REGIONWIDTH] == size(shiftedSubRegion_01))
        REGION_MATRIX_01(:, :, n) = shiftedSubRegion_01;
    else
        paddedRegion = zeros(REGIONHEIGHT, REGIONWIDTH);
        paddedRegion( 1 + max(0, 1-ymin1(n)) : REGIONHEIGHT - max(0, ymax1(n)-imageHeight),...
            1 + max(0, 1 - xmin1(n)) : REGIONWIDTH - max(0, xmax1(n) - imageWidth)) = shiftedSubRegion_01;
        REGION_MATRIX_01(:, :, n) = paddedRegion;
    end

    % If the dimensions of the second shifted region are less than the specified
    % dimensions of the interrogation regions, then Pad the image subregion
    % with zeros so that its final dimensions are consistent with those
    % specified.
    if all([REGIONHEIGHT, REGIONWIDTH] == size(shiftedSubRegion_02))
        REGION_MATRIX_02(:, :, n) = shiftedSubRegion_02;
    else
        paddedRegion = zeros(REGIONHEIGHT, REGIONWIDTH);
        paddedRegion( 1 + max(0, 1 - ymin2(n)) : REGIONHEIGHT - max(0, ymax2(n) - imageHeight),...
            1 + max(0, 1 - xmin2(n)) : REGIONWIDTH - max(0, xmax2(n) - imageWidth)) = shiftedSubRegion_02;
        REGION_MATRIX_02(:, :, n) = paddedRegion;
    end
    
end

end






