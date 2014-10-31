function rpcErrorAnalysisMonteCarlo(JOBFILE, IMDIR, SAVEPATH, IMAGEFILEPATH, NDIGITS, PARAMETERSPATH, SPATIALWINDOW, IMAGESPECTRALFILTER, CONSTANTS, SAVEDATA, NPROCESSORS)

% Read the job file (this is just for saving at the end of the function)
JobFile = JOBFILE;

% Default to not saving data
if nargin < 11
    SAVEDATA = 0;
end

% Specify numbering format based on number of digits
numberFormat = ['%0' num2str(NDIGITS) '.0f'];

% try
    
% Load in known transformation parameters (used for error analysis)
load(PARAMETERSPATH);
trueRotation = Parameters.Rotation;
trueScaling = Parameters.Scaling;
concentration = Parameters.Concentration;
trueTransformationMatrix = Parameters.Tforms;
trueTranslationX = Parameters.TranslationX;
trueTranslationY = Parameters.TranslationY;
shearX = Parameters.ShearX;
shearY = Parameters.ShearY;

% Noise mean and standard deviation fractions
noiseMean = CONSTANTS.NoiseMean;
noiseStd = CONSTANTS.NoiseStd;

% Number of images
nImages = length(trueTranslationX); 

% Start and stop files
startFile = 1;
stopFile = nImages;

% Size of raw images
imageHeight = Parameters.ImageHeight;
imageWidth = Parameters.ImageWidth;

% Initialize vectors to hold estimates
estimatedTranslationX = zeros(nImages, 1);
estimatedTranslationY = zeros(nImages, 1);

% Initialize vectors to hold absolute translation errors
translationAbsErrorX = zeros(nImages, 1);
translationAbsErrorY = zeros(nImages, 1);

% Initialize vector to hold peak ratios
peakRatio = zeros(nImages, 1);

% Load the image matrix.
imageData = load(IMAGEFILEPATH);
imageMatrix1 = double(imageData.imageMatrix1);
imageMatrix2 = double(imageData.imageMatrix2);

% Max value of the images
maxVal = double(intmax(class(imageData.imageMatrix1)));

% Make noise matrices. The 2.8 corresponds to the multiple of the standard
% deviation corresponding to a 99.5% coverage factor.
noiseMatrix1 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(imageMatrix1));
noiseMatrix2 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(imageMatrix2));

% Do the processing
if NPROCESSORS > 1

    % Loop through images
    parfor k = 1: nImages

    % Read in the first image in the pair 
    IMAGE1 = imageMatrix1(:, :, k) + noiseMatrix1(:, :, k);

    % Read in the second image in the pair
    IMAGE2 = imageMatrix2(:, :, k) + noiseMatrix2(:, :, k);

    % Window the images
        ImgWindowed1 = SPATIALWINDOW .* IMAGE1;           
        ImgWindowed2 = SPATIALWINDOW .* IMAGE2;

    % Calculate RPC between images. The 0 at the end means "Only take the primary peak.     
        [estimatedTranslationY(k), estimatedTranslationX(k), correlationPlane] = RPC(ImgWindowed1, ImgWindowed2, IMAGESPECTRALFILTER, 0);


    end % End (parfor k = 1 : nImages )

else % Else if nProcessors == 1

    for k = 1: nImages

        % Read in the first image in the pair 
        IMAGE1 = imageMatrix1(:, :, k) + noiseMatrix1(:, :, k);

        % Read in the second image in the pair
        IMAGE2 = imageMatrix2(:, :, k) + noiseMatrix2(:, :, k);

        % Window the images
        ImgWindowed1 = SPATIALWINDOW .* IMAGE1;           
        ImgWindowed2 = SPATIALWINDOW .* IMAGE2;            

        % Calculate RPC between images. The 0 at the end means "only take the primary peak."     
        [estimatedTranslationY(k), estimatedTranslationX(k), correlationPlane] = RPC(ImgWindowed1, ImgWindowed2, IMAGESPECTRALFILTER, 0);


    end % End (for k = 1 : nImages)

end % End ( if nProcessors > 1)
    
    
% Calculate translation errors (i.e. difference between the calculated and true values for translation
translationErrorX =  trueTranslationX  - estimatedTranslationX;
translationErrorY =  trueTranslationY - estimatedTranslationY;
    
% Extract constants from structure
spatialWindowFraction = CONSTANTS.SpatialWindowFraction;
spatialRPCdiameter = CONSTANTS.SpatialRPCdiameter;

% rename image directory variable
imageDirectory = IMDIR;

% Save outputs to disk
save( SAVEPATH, 'JobFile', 'concentration', 'trueRotation', ...
'trueScaling', 'trueTranslationX', 'trueTranslationY', ...
'estimatedTranslationX', 'estimatedTranslationY', ...
 'translationErrorX', 'translationErrorY', 'peakRatio', 'spatialWindowFraction', ...
'spatialRPCdiameter', 'imageHeight', 'imageWidth', 'shearX', 'shearY');




end




