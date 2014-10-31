
function [y, x, v, u, r, s, uVal, vVal, rVal] = addNoise_fullField(FILEPATHS, JOBFILE)

% tic
% profile on;

% Rename variables
JobFile = JOBFILE;
FilePaths = FILEPATHS;

% Number of PIV passes
% Take the minimum of the requested number of passes and the number of pass
% parameter structures specified.
numberOfPasses = min(JobFile.JobOptions.NumberOfPasses, length(JobFile.Parameters.Processing));

% Simulate beam?
simulateBeam = JobFile.JobOptions.SimulateBeam;

% Simulate noise? 
simulateNoise = JobFile.JobOptions.SimulateNoise;

% Angle by which to rotate images.
if isempty(JOBFILE.JobOptions.ImageRotationAngle)
    rotationAngle = 0;
else
    rotationAngle = JOBFILE.JobOptions.ImageRotationAngle;
end

% Read in the clean images.
image1clean = imread(FilePaths.FirstImagePath);
image2clean = imread(FilePaths.SecondImagePath);

% Image dimensions
[imageHeight, imageWidth] = size(image1clean);
imageSize = [imageHeight, imageWidth];

% Beam profile standard deviation, so that FWHM width
% is equal to half the image width.
if simulateBeam
    beamStd = imageWidth / (2 * sqrt(2 * log(2)));
    beamCenter = imageWidth / 2;
    imageGridX = 1 : imageWidth;
    beamProfile = exp(-(imageGridX - beamCenter).^2 / (2 * beamStd.^2));
    beamImage = repmat(beamProfile, [imageHeight, 1]);
else
    beamImage = ones(imageHeight, imageWidth);
end

% Attenuate the particle intensities according to the specified
% beam profile.
image1Beam = double(image1clean) .* beamImage;
image2Beam = double(image2clean) .* beamImage;

% Max value of the images
maxVal = double(intmax(class(image1clean)));

% Noise mean and standard deviation fractions
noiseMean = JobFile.Parameters.Noise.Mean;
noiseStd = JobFile.Parameters.Noise.Std;

% Make noise matrices. The 2.8 corresponds to the multiple of the standard
% deviation corresponding to a 99.5% coverage factor.
% There is really no reason to be adding noise to the images;
% this should be done in image generation. But I already generated
% a lot of images and don't want to re-gen them every time
% I change noise parameters.
if simulateNoise
    noiseMatrix1 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(image1Beam));
    noiseMatrix2 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(image2Beam));
    % Background intensity level
    noiseBackground = JobFile.Parameters.Noise.Background;
else
    noiseMatrix1 = zeros(imageHeight, imageWidth);
    noiseMatrix2 = zeros(imageHeight, imageWidth);
    noiseBackground = 0;
end

% Add noise to the images. Rotate the images according to the specifed
% angle, if requested.
image1 = imrotate(double(image1Beam), rotationAngle, 'bicubic') + noiseMatrix1 + (noiseBackground * maxVal) * beamImage;
image2 = imrotate(double(image2Beam), rotationAngle, 'bicubic') + noiseMatrix2 + (noiseBackground * maxVal) * beamImage;

% Saturate and zero pixels
image1(image1 > maxVal) = maxVal;
image2(image2 > maxVal) = maxVal;
image1(image1 < 0) = 0;
image2(image2 < 0) = 0;

% Loop over the passes.
for p = 1 : numberOfPasses;

% Interrogation Region Parameters
regionHeight = JobFile.Parameters.Processing(p).InterrogationRegion.Height;
regionWidth = JobFile.Parameters.Processing(p).InterrogationRegion.Width;     
    
% Grid parameters
gridSpacingX = JobFile.Parameters.Processing(p).Grid.Spacing.X;
gridSpacingY = JobFile.Parameters.Processing(p).Grid.Spacing.Y;

% Make sure the grid buffer is at least half the size of the interrogation region
gridBufferY = max(JobFile.Parameters.Processing(p).Grid.Buffer.Y, ceil(regionHeight / 2));
gridBufferX = max(JobFile.Parameters.Processing(p).Grid.Buffer.X, ceil(regionWidth / 2));

% Generate the list of coordinates that specifies the (X, Y) centers of all of the interrogation regions 
[ gx_01, gy_01 ] = gridImage(imageSize, [gridSpacingY gridSpacingX], gridBufferY, gridBufferX);

% If the pass number is greater than one, i.e., if at least one pass has
% finished, then interpolate the velocity field from the previous pass
% onto the grid for the current pass.
% Round the grid shift values so that grid points are shifted
% from integer coordinates to integercoordinates
if p > 1
    gridShiftX = round(interp2(X{p-1}, Y{p-1}, UVAL{p-1}, gx_01, gy_01, 'linear', 0));
    gridShiftY = round(interp2(X{p-1}, Y{p-1}, VVAL{p-1}, gx_01, gy_01, 'linear', 0));
else
    gridShiftX = zeros(size(gx_01));
    gridShiftY = zeros(size(gx_01));
end

% Grid the second image, taking into account the discrete window offsets.
[ gx_02, gy_02 ] = gridImage(imageSize, [gridSpacingY gridSpacingX], gridBufferY, gridBufferX, gridShiftY, gridShiftX);

% Determine the size of the grid (number of grid points)..
[numRows, numColumns] = size(gx_01);

% Reshape the grid point arrays for the first image into vectors
gridPointsX_01 = gx_01(:);
gridPointsY_01 = gy_01(:);

% Reshape the grid point arrays for the second image into vectors
gridPointsX_02 = gx_02(:);
gridPointsY_02 = gy_02(:);

% Determine the number of interrogation regions to be correlated
nRegions = length(gridPointsX_01);

% Apodization window parameters.
spatialWindowFraction = JobFile.Parameters.Processing(p).InterrogationRegion.SpatialWindowFraction;
fmiWindowSize = JobFile.Parameters.Processing(p).InterrogationRegion.FMIWindowSize;
fmiWindowType = JobFile.Parameters.Processing(p).InterrogationRegion.FMIWindowType;

% FFT Parameters
fftSize = JOBFILE.Parameters.Processing(p).FFTSize;
fftHeight = fftSize(1);
fftWidth = fftSize(2);

% Image resampling parameters
numberOfRings = JobFile.Parameters.Processing(p).Resampling.NumberOfRings;
numberOfWedges = JobFile.Parameters.Processing(p).Resampling.NumberOfWedges;
rMin = JobFile.Parameters.Processing(p).Resampling.MinimumRadius;
rMax = min(fftHeight, fftWidth) / 2 - 1;

% Correlation parameters
correlationMethod = JobFile.Parameters.Processing(p).Correlation.Method;
spatialRPCDiameter = JobFile.Parameters.Processing(p).Correlation.SpatialRPCDiameter;
fmiRpcDiameter = JobFile.Parameters.Processing(p).Correlation.FMCDiameter; 

% Create the gaussian intensity window to be applied to the the raw image interrogation regions
% spatialWindow = gaussianWindowFilter( [regionHeight regionWidth], spatialWindowFraction, 'fraction');
spatialWindow = gaussianWindowFilter_old( [regionHeight regionWidth], spatialWindowFraction, 'fraction');

% Determine the FMI window type.
isHann1 = ~isempty(regexpi(fmiWindowType, 'hann1'));
isHann2 = ~isempty(regexpi(fmiWindowType, 'hann2'));
isGaussianSkew = ~isempty(regexpi(fmiWindowType, 'gauss_skew'));

% Create the FMI Window
if isHann1
    fmiWindow1D = hann1(numberOfRings, [fmiWindowSize(1) fmiWindowSize(2)], fmiWindowSize(3));
    fmiWindow = repmat(fmiWindow1D, numberOfWedges, 1);
elseif isHann2
    fmiWindow = hann2([numberOfWedges, numberOfRings], fmiWindowSize(1));
elseif isGaussianSkew
    fmiWindow1D = gaussianWindowFilter_asymmetric(numberOfRings, fmiWindowFraction);
    fmiWindow = repmat(fmiWindow1D, numberOfWedges, 1);
else
%     fmiWindow = gaussianWindowFilter([numberOfWedges, numberOfRings], fmiWindowSize, 'fraction');
    fmiWindow = gaussianWindowFilter_old([numberOfWedges, numberOfRings], fmiWindowSize, 'fraction');
end

% Create the gaussian spectral energy filter be applied to the raw image correlation
imageSpectralFilter = spectralEnergyFilter(regionHeight, regionWidth, spatialRPCDiameter); 

% Create the FMI spectral filter (i.e. the FMI RPC filter).
fmiSpectralFilter = spectralEnergyFilter(numberOfWedges, numberOfRings, fmiRpcDiameter); 

% Preallocate memory for the vectors to hold the estimates of translation, rotation, and scaling.
estimatedTranslationY = zeros(nRegions, 1);
estimatedTranslationX = zeros(nRegions, 1);
estimatedRotation = zeros(nRegions, 1); 
estimatedScaling = zeros(nRegions, 1);
fmiTranslationY = zeros(nRegions, 1);
fmiTranslationX = zeros(nRegions, 1);

% Make a matrix of the subregion coordinates.
% Do this only once to increase speed (meshgrid is slow).
[xImage, yImage] = meshgrid(1 : regionWidth, 1 : regionHeight);

% FFT Spectrum coordinates
[xSpectrum, ySpectrum] = meshgrid(1:fftWidth, 1:fftHeight);

% Figure out the log polar resampling coordinates.
[xLP, yLP] = LogPolarCoordinates([fftHeight, fftWidth], numberOfWedges, numberOfRings, rMin, rMax, 2 * pi);

% ftFilter = fftshift(makeFFTFilter([imageHeight, imageWidth], 1));
ftFilter = ones([numberOfWedges, numberOfRings]);

% Initialize FMC peak height ratio vector.
fmcPeakRatio = zeros(nRegions, 1);

% Initialize RPC peak height ratio vector.
spatialPeakRatio = zeros(nRegions, 1);

% Initialize affine parameters matrix
affineParameters = zeros(nRegions, 6);

% Determine which correlation type to use. Fmc or RPC
isFmc = strcmpi(correlationMethod, 'fmc');
isRpc = strcmpi(correlationMethod, 'rpc');
isScc = strcmpi(correlationMethod, 'scc');

% Determine the leftmost and rightmost columns of the interrogation regions
% (first image)
xMin_01 = gridPointsX_01 - ceil( regionWidth / 2 ) + 1;
xMax_01 = gridPointsX_01 + floor( regionWidth / 2 );

% Determine the top and bottom rows of the interrogation regions
% (first image)
yMin_01 = gridPointsY_01 - ceil( regionHeight / 2 ) + 1;
yMax_01 = gridPointsY_01 + floor ( regionWidth / 2 );

% Determine the leftmost and rightmost columns of the interrogation regions
% (second image)
xMin_02 = gridPointsX_02 - ceil( regionWidth / 2 ) + 1;
xMax_02 = gridPointsX_02 + floor( regionWidth / 2 );

% Determine the top and bottom rows of the interrogation regions
% (second image)
yMin_02 = gridPointsY_02 - ceil( regionHeight / 2 ) + 1;
yMax_02 = gridPointsY_02 + floor ( regionWidth / 2 );

% Initialze a matrix for holding a stack of subregions.
regionMatrix1 = zeros(regionHeight, regionWidth, nRegions);
regionMatrix2 = zeros(regionHeight, regionWidth, nRegions);

% Create stacks of the subregions.
% This is done so that individual subregions (rather than entire images)
% can be passed to the parfor loop. 
for k = 1 : nRegions
    regionMatrix1(:, :, k) = image1( max([1, yMin_01(k)]) : min([imageHeight, yMax_01(k)]) , max([1, xMin_01(k)]) : min([imageWidth, xMax_01(k)]) );
    regionMatrix2(:, :, k) = image2( max([1, yMin_02(k)]) : min([imageHeight, yMax_02(k)]) , max([1, xMin_02(k)]) : min([imageWidth, xMax_02(k)]) );
end

% COMPILED = ~ismac;
COMPILED = 0;

% Start a timer
t = tic;
parfor k = 1 : nRegions
% for k = 417 : 417
    
%     disp([num2str(k) ' of ' num2str(nRegions)]);
% for k = 104 : 104;
    % Extract the subregions from the subregion stacks.
    subRegion1 = regionMatrix1(:, :, k);
    subRegion2 = regionMatrix2(:, :, k);
   
%     % Apply the spatial window to the subregions
%     windowedRegion1 = spatialWindow .* subRegion1;
%     windowedRegion2 = spatialWindow .* subRegion2;

    % Perform FMC processing. 
    if isFmc
        % Perform the FMC correlation.
        [estimatedTranslationY(k), estimatedTranslationX(k),...
        estimatedRotation(k), estimatedScaling(k), ...
        fmcPeakRatio(k), spatialPeakRatio(k)] = ...
        ...
        FMC(subRegion1, subRegion2, spatialWindow, imageSpectralFilter,...
        fmiWindow, fmiSpectralFilter, ...
        xImage, yImage, ...
        xSpectrum, ySpectrum,...
        xLP, yLP, rMin, rMax, COMPILED);
    
    % Perform RPC analysis 
    % The zero in this input means "Do not search multiple peaks,"
    % i.e., use only the primary peak.
    elseif isRpc
        [estimatedTranslationY(k), estimatedTranslationX(k), rpcPlane]...
            = RPC(spatialWindow .* subRegion1, spatialWindow .* subRegion2,...
            imageSpectralFilter); 
        
        % Measure the peak height ratio
        if COMPILED
            spatialPeakRatio(k) = measurePeakHeightRatio_mex(rpcPlane);
        else
            spatialPeakRatio(k) = measurePeakHeightRatio(rpcPlane);
        end
        
    
    % Perform SCC analysis.
    elseif isScc
        [estimatedTranslationY(k), estimatedTranslationX(k), spatialPeakRatio(k)]...
            = SCC(spatialWindow .* subRegion1, spatialWindow .* subRegion2);
    end    
end % end for k = 1 : nPairs

disp(['Correlation times (pass ' num2str(p) '): ' num2str(toc(t)) ' sec' ]);

% profile viewer
 % Rearrange the grid matrices
x = (reshape(gridPointsX_01, numRows, numColumns));
y = (reshape(gridPointsY_01, numRows, numColumns));

% Reshape the results matrices
% Add the applied DWO grid shift to the translation matrices, so that the output
% translations aren't just fractions of pixels!!
TRANSLATIONX = (reshape(estimatedTranslationX, numRows, numColumns)) + gridShiftX;
TRANSLATIONY = (reshape(estimatedTranslationY, numRows, numColumns)) + gridShiftY;
ROTATION = reshape(estimatedRotation, numRows, numColumns);
SCALING = reshape(estimatedScaling, numRows, numColumns);
spatialPeakRatio = reshape(spatialPeakRatio, numRows, numColumns);

%  Now flip everything back across the X-axis.
x = flipud(x);
y = flipud(y);
u = flipud(TRANSLATIONX);
v = flipud(TRANSLATIONY);
r = flipud(ROTATION);
s = flipud(SCALING);

% Calculate gradient-based vorticity
[uVal, vVal, rVal] = validateField_2d3c(x, y, u, v, r);

% Save the results from each pass to the output variable structures.
X{p} = x;
Y{p} = y;
U{p} = u;
V{p} = v;
R{p} = r;
S{p} = s;
UVAL{p} = uVal;
VVAL{p} = vVal;
RVAL{p} = rVal;
FMC_PEAK_RATIO{p} = fmcPeakRatio;
SPATIAL_PEAK_RATIO{p} = spatialPeakRatio;

end % End for p = 1 : numberOfPasses

% Save the results
save(FilePaths.OutputFilePath, ...
    'X', 'Y', 'U', 'V', 'R', 'S',...
    'UVAL', 'VVAL', 'RVAL', ...
    'FMC_PEAK_RATIO', 'SPATIAL_PEAK_RATIO',...
    'FilePaths', 'JobFile');

toc(t);

end






















