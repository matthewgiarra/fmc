function [ANGLEOUT, SCALINGOUT, TRANSLATIONX, TRANSLATIONY, RPC_PEAK_HEIGHT_RATIO, RPC_PEAK_HEIGHT, RPC_PEAK_DIAMETER, PEAKUSED] = resolveRotationAmbiguity(IMAGE_01, IMAGE_02, SPATIAL_WINDOW, SPECTRAL_FILTER, ANGLE_IN, SCALING_IN, Y_IMAGE, X_IMAGE, MULTIPEAK, DIFFERENCE_METHOD, COMPILED)
% [ANGLEOUT, SCALINGOUT, TRANSLATIONX, TRANSLATIONY, RPC_PEAK_HEIGHT_RATIO, RPC_PEAK_HEIGHT, RPC_PEAK_DIAMETER, PEAKUSED] = resolveRotationAmbiguity(IMAGE_01, IMAGE_02, SPATIAL_WINDOW, SPECTRAL_FILTER, ANGLE_IN, SCALING_IN, Y_IMAGE, X_IMAGE, MULTIPEAK, DIFFERENCE_METHOD, COMPILED)
% This function aligns the images by the calculated rotation angle and
% finds the translation relating them.
%
% INPUTS:
%   IMAGE_01 = [M x N] grayscale image specifying first interrogation 
%       region (raw, un-rotated, un-windowed). M and N are the number of 
%       rows and columns in the image (in pixels).
%
%   IMAGE_02 = [M x N] grayscale image specifying first interrogation
%       region (raw, un-rotated, un-windowed).
%
%   SPATIAL_WINDOW = [M x N] matrix specifying apodization window applied 
%       to IMAGE_01 and IMAGE_02 (e.g., Gaussian window).
%   
%   SPECTRAL_FILTER = [M x N] matrix specifying the spectral filter applied
%       to the spatial correlation plane relating the rotated/scaled copies of
%       IMAGE_01 and IMAGE_02.
%
%   ANGLEIN = Initial estimate of rotation angle (in radians) relating
%       IMAGE_01 and IMAGE_02, usually calculated using the function FMIRPC.m
%
%   SCALING_IN = Initial estimate of scaling factor(unitless) relating
%       IMAGE_01 and IMAGE_02, usually calculated using the function FMIRPC.m
%
%   Y_IMAGE = [M x N] matrix specifying the row coordinates of IMAGE_01 and
%       IMAGE_02, shifted such that the center of rotation
%       (i.e., (0,0) coordinate) is coincident with the geometric centroid of
%       the image. This is passed as an argument rather than calculated
%       inside the function because MESHGRID (and similar functions) are
%       rather expensive and really need not be calculated every time this
%       funciton is run.
%
%   X_IMAGE = [M x N] matrix specifying the column coordinates of IMAGE_01 and
%       IMAGE_02, shifted such that the center of rotation
%       (i.e., (0,0) coordinate) is coincident with the geometric centroid of
%       the image. 
%
%   MULTIPEAK = Boolean flag specifying whether or not to check the
%       rotations/scalings that corresponded to multiple peaks in the FMC
%       plane. Note that in order for this flag to have any effect on this
%       function, the same flag must have been enabled in the function
%       FMIRPC.m, so that multiple rotation angles / scaling factors were
%       output.
%
%   DIFFERENCE_METHOD = Differencing method for rotating / scaling image
%       pair (central difference, forward difference, or backward difference).
%       This option is specifed in the following way:
%
%       DIFFERENCE_METHOD = 1 (default) specifies central difference, wherein
%       IMAGE_01 is rotated by (ANGLE_IN/2) and scaled by
%       (SCALING_IN)^(1/2), and IMAGE_02 is rotated by -(ANGLE_IN/2) and scaled by
%       (SCALING_IN)^(-1/2).
%
%       DIFFERENCE_METHOD = 2 specifies forward difference, wherein
%       IMAGE_01 is rotated by ANGLE_IN and scaled by
%       (SCALING_IN), and IMAGE_02 is left untouched (i.e., neither rotated
%       nor scaled). 
%
%       DIFFERENCE_METHOD = 3 specifies backward difference, wherein
%       IMAGE_01 is left untouched (i.e., neither rotated
%       nor scaled), and IMAGE_02 is rotated by -(ANGLE_IN/2) and scaled by
%       (SCALING_IN)^(-1).
%
%   COMPILED = Boolean flag that specifies whether or not to use compiled 
%       interpolation codes to perform the rotation / scaling of images.
%
% OUTPUTS:
%   ANGLEOUT = Rotation angle (radians) resulting in highest peak 
%   in the phase correlation between the rotated and scaled image pair.
%   
%   SCALINGOUT = Scaling factor (unitless) resulting in highest peak
%       in the phase correlation between the rotated and scaled image pair.
%   
%   TRANSLATIONX = RPC estimate of the horizontal (column)
%       translation (in pixels; non-integer) relating the pair of images
%       that were rotated by ANGLEOUT and scaled by SCALINGOUT.
%   
%   TRANSLATIONY = RPC estimate of the vertical (row)
%       translation (in pixels; non-integer) relating the pair of images
%       that were rotated by ANGLEOUT and scaled by SCALINGOUT. Note that
%       these translations assume a left-handed coordinate system with the
%       positive-vertical axis pointing from the top to the bottom of the
%       image (e.g., positive y-translation means the pattern moved
%       downward).
%
%   RPC_PEAK_HEIGHT_RATIO = Peak height ratio of the RPC plane, calculated
%       once for each FMC peak used. 
%
%   RPC_PEAK_HEIGHT = Height of the peak of the RPC plane, calculated once
%   for each FMC peak used.
%
%   RPC_PEAK_DIAMETER = Diameter of the RPC peak in pixels. See the
%   function subpixel.m for details on this calculation.
%
%   PEAKUSED = Number of the FMC peak resulting in the greatest value of
%       RPC_PEAK_HEIGHT.
%  
% SEE ALSO:
%   FMIRPC, subpixel


% Default to not using compiled codes.
if nargin < 11
    COMPILED = 0;
end;

% Default to central difference method
% DIFFERENCE_METHOD = 1 is central difference
% DIFFERENCE_METHOD = 2 is forward difference
% DIFFERENCE_METHOD = 3 is backward difference
if nargin < 10
    DIFFERENCE_METHOD = 1;
end

% Default to not using multiple peaks
if nargin < 9
    MULTIPEAK = 0;
end

% Number of FMC peaks that were used to estimate rotation and scaling.
nPeaks = 1 + (length(ANGLE_IN) - 2) * MULTIPEAK;

% Initialize the peak ratio vector. 
peakRatio = zeros(nPeaks, 1);

% Initialize the peak height vector
peakHeight = zeros(nPeaks, 1);

% Initialize transformation parameters.
angleOut = zeros(nPeaks, 1);
tx = zeros(nPeaks, 1);
ty = zeros(nPeaks, 1);

% % Loop over the peak number (p) 
for p = 1 : nPeaks;     
    [angleOut(p), tx(p), ty(p), peakHeight(p), peakRatio(p), peakDiameter(p)] = ...
    testAffineParameters(double(IMAGE_01), double(IMAGE_02), double(SPATIAL_WINDOW), double(SPECTRAL_FILTER), ANGLE_IN(p), SCALING_IN(p), X_IMAGE, Y_IMAGE, DIFFERENCE_METHOD, COMPILED);
end

% Determine which of the rotations/scalings corresponding to the n FMC 
% peaks had the highest RPC correlation peak.
[~, PEAKUSED] = max(peakHeight);

% % Save variables to output.
ANGLEOUT = angleOut(PEAKUSED);
SCALINGOUT = SCALING_IN(PEAKUSED);
TRANSLATIONX = tx(PEAKUSED);
TRANSLATIONY = ty(PEAKUSED);
RPC_PEAK_HEIGHT_RATIO = peakRatio(PEAKUSED);
RPC_PEAK_HEIGHT = peakHeight(PEAKUSED);
RPC_PEAK_DIAMETER = peakDiameter(PEAKUSED);

end


function [angleOut, tx, ty, peakHeight, peakRatio, peakDiameter] = testAffineParameters(image1, image2, spatialWindow, spectralFilter, angleIn, scalingIn, xImage, yImage, differenceMethod, COMPILED)

% Invert the rotation angle
th = -angleIn;

% Similarity matrices
[matrix_11, matrix_21] = makeSimilarityMatrices(th, scalingIn, differenceMethod);

% These lines transform the images. The transform function is written so
% that if an identity matrix is input, it quickly outputs the original
% image, so it shouldn't be inefficient to do this outside the "if"
% statements.
% This transforms the first image.
imageTest_1 = transformImage(image1, xImage, yImage, matrix_11, COMPILED);

% This transforms the second image.
imageTest_2 = transformImage(image2, xImage, yImage, matrix_21, COMPILED);

% This is does the RPC correlation between the aligned first and second
% images that were aligned using the originally-calculated rotation angle.
[ty1, tx1, spatialCorr, peakHeight1, peakDiameter1] = RPC(spatialWindow .* imageTest_1, spatialWindow .* imageTest_2, spectralFilter, COMPILED);

% This calculates the peak ratio and peak height for the correlation
% between the "central difference" correlation for the
% originally-calculated rotation angle.
[peakRatio1, ~] = measurePeakHeightRatio(spatialCorr, COMPILED);
% peakRatio1 = calculate_mutual_information(imageTest_1, imageTest_2);

% Check for the other angle. Keep the original variable names to save speed
% on memory allocation.
 
% Similarity matrices
[matrix_12, matrix_22] = makeSimilarityMatrices(th + pi, scalingIn, differenceMethod);

% This is the FIRST image rotated by half the originally calculated
% rotation angle + pi and scaled by the square root of the originally calculated
% scaling factor, all in the FORWARD difference direction. 
imageTest_1 = transformImage(image1, xImage, yImage, matrix_12, COMPILED);

% This is the SECOND image rotated by half the originally calculated
% rotation angle + pi and scaled by the square root of the originally calculated
% scaling factor, all in the BACKWARD difference direction. 
imageTest_2 = transformImage(image2, xImage, yImage, matrix_22, COMPILED);

% This is does the RPC correlation between the forward-transformed first
% image and the backward-transformed second image, for the originally
% calculated rotation angle + pi.
[ty2, tx2, spatialCorr, peakHeight2, peakDiameter2] = RPC(spatialWindow .* imageTest_1, spatialWindow .* imageTest_2, spectralFilter, COMPILED);

% This calculates the peak ratio and peak height for the correlation
% between the "central difference" correlation for the
% originally-calculated rotation angle + pi.
[peakRatio2, ~] = measurePeakHeightRatio(spatialCorr, COMPILED);
% peakRatio2 = calculate_mutual_information(imageTest_1, imageTest_2);

% If backward difference was used, then the vector needs to be rotated and
% scaled to transform it into forward-time velocity.
if differenceMethod == 3
    % Set up some temporary variables so we don't overwrite the important
    % stuff.
    tx1_temp = tx1;
    ty1_temp = ty1;
    tx2_temp = tx2;
    ty2_temp = ty2;

    % This corrects the measured translation to account for 
    % having rotated the second image by the originally calculated
    % rotation angle.
    tx1 = s * (tx1_temp * cos(th) - ty1_temp * sin(th));
    ty1 = s * (tx1_temp * sin(th) + ty1_temp * cos(th)) ;

    % This corrects the measured translation to account for 
    % having rotated the second image by the originally calculated
    % rotation angle + pi.
    tx2 = s * (tx2_temp * cos(th + pi) - ty2_temp * sin(th + pi));
    ty2 = s * (tx2_temp * sin(th + pi) + ty2_temp * cos(th + pi));
end

% Decide between rotation angles
% if peakRatio1 >= peakRatio2 % Originally calculated angle was better
if peakHeight1 >= peakHeight2
    rotationAngle =  -th;
    tx = tx1;
    ty = ty1;
    peakRatio = peakRatio1;
    peakHeight = peakHeight1;
    peakDiameter = peakDiameter1;
    
else % Calculated angle + pi was better
    rotationAngle = -th + pi;
    tx = tx2;
    ty = ty2;
    peakRatio = peakRatio2;
    peakHeight = peakHeight2;
    peakDiameter = peakDiameter2;
end

% Make the rotation interval -pi to pi
% and save resultant angle to output variable. 
angleOut = normalizeAngle(rotationAngle, 0);

end




