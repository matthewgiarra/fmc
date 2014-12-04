function [DISPARITY_X, DISPARITY_Y] = calculateTransformedRegionDisparity(REGION_01, REGION_02, TRANSLATION_Y, TRANSLATION_X, ROTATION_ANGLE, SCALING_FACTOR, YGRID, XGRID, CONFIDENCE_INTERVAL, DIFFERENCE_METHOD, COMPILED)

% Determine height and width of interrogation regions.
[region_height, region_width] = size(REGION_01);

% % Calculate the horizontal and vertical coordinate of the
% geometric centroid of the region.
xc = region_width / 2 - 0.5;
yc = region_height / 2 - 0.5;

% Shift coordinates to put origin at region center.
x_grid = (XGRID - 1) - xc;
y_grid = (YGRID - 1) - yc;

% Calculate similarity matrices for the rotation and translation
[matrix_01, matrix_02] = makeSimilarityMatrices_homogeneous(ROTATION_ANGLE, SCALING_FACTOR, TRANSLATION_Y, TRANSLATION_X, DIFFERENCE_METHOD);

% Transform region 1
region_01_transformed = transformImage_homogeneous(REGION_01, x_grid, y_grid, matrix_01, COMPILED); 

% Transform region 2
region_02_transformed = transformImage_homogeneous(REGION_02, x_grid, y_grid, matrix_02, COMPILED); 

% Calculate the disparity between the images.
[DISPARITY_X, DISPARITY_Y] = regionDisparity(region_01_transformed, region_02_transformed, CONFIDENCE_INTERVAL);

end