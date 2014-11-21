function fmiErrorAnalysisMonteCarlo(JOBFILE, SAVEPATH, IMAGEFILEPATH, PARAMETERSPATH, SPATIALWINDOW, IMAGESPECTRALFILTER, CONSTANTS, NPROCESSORS)

% Rename the Jobfile variable.
JobFile = JOBFILE;

% Flag specifying whether to run compiled codes.
run_compiled = JobFile.JobOptions.RunCompiled;

% Load the image matrix.
imageData = load(IMAGEFILEPATH);

% Matrices containing the image data for the first and second image in each pair.
imageMatrix1 = double(imageData.imageMatrix1);
imageMatrix2 = double(imageData.imageMatrix2);

% Number of images
nImages = size(imageMatrix1, 3);

% Window fraction for FMI
fmiWindowFraction = CONSTANTS.FmiWindowFraction;
fmiRPCDiameter = CONSTANTS.FmiRPCdiameter;

% Minimum and maximum number of wedges (angular sampling resolution) in log polar resampling.
NWmin = CONSTANTS.NWmin;
NWmax = CONSTANTS.NWmax;

% Minimum and maximum number of rings (radial sampling resolution) in log polar resampling.
NRmin = CONSTANTS.NRmin;
NRmax = CONSTANTS.NRmax;

% FMIWindow type
fmiWindowType = CONSTANTS.FMIWindowType;

% Dimensions of image FFTs
fftSize = CONSTANTS.FFTSize;

% Load in known transformation parameters (used for error analysis)
load(PARAMETERSPATH);
trueRotation = Parameters.Rotation;
trueScaling = Parameters.Scaling;
trueTranslationX = Parameters.TranslationX;
trueTranslationY = Parameters.TranslationY;

% Number of wedges in the FM- image
numWedgesOdd = ceil( (NWmin + (NWmax - NWmin) * rand(nImages, 1)));
numWedges = numWedgesOdd + mod(numWedgesOdd, 2); % Force to even integer

% Number of rings in the FM-image
numRingsOdd = ceil( (NRmin + (NRmax - NRmin) * rand(nImages, 1)));
numRings = numRingsOdd + mod(numRingsOdd, 2); % Force to even integer

% Size of raw images
imageHeight = Parameters.ImageHeight;
imageWidth = Parameters.ImageWidth;

% Height and width of the image that's going to be resampled.
% This will be the height and width of the FFT if doing FMC
spectrum_height = fftSize(1);
spectrum_width  = fftSize(2);

% Minimum and maximum radii of log-polar resampling of the FFT magnitude. 
rMax = min(spectrum_height, spectrum_width) / 2 - 1;
rMin = CONSTANTS.RMin;

% Make the image coordinates
[xImage, yImage] = meshgrid(1:imageWidth, 1:imageHeight);

% Figure out the log polar resampling coordinates.
[xLP, yLP] = LogPolarCoordinates([spectrum_height, spectrum_width], numWedges(1), numRings(1), rMin, rMax, 2*pi);

% Initialize vectors to hold estimates
estimatedRotation = zeros(nImages, 1);
estimatedScaling = zeros(nImages, 1);
estimatedTranslationX = zeros(nImages, 1);
estimatedTranslationY = zeros(nImages, 1);

% If there was no distribution of num wedges specified then just make a
% single FMI window and spectral filter
uniformWindowDistribution = min(numWedges) == max(numWedges) && min(numRings) == max(numRings) ;

% Determine the window type.
isHann1 = ~isempty(regexpi(fmiWindowType, 'hann1'));
isHann2 = ~isempty(regexpi(fmiWindowType, 'hann2'));
isGaussianSkew = ~isempty(regexpi(fmiWindowType, 'gauss_skew'));

if uniformWindowDistribution
    % Create the gaussian windows
    if isHann1
        fmiWindowUniform1D = hann1(numRings(1), [fmiWindowFraction(1) fmiWindowFraction(2)], fmiWindowFraction(3));
        fmiWindowUniform = repmat(fmiWindowUniform1D, numWedges(1), 1);
    elseif isHann2
        fmiWindowUniform = hann2([numWedges(1), numRings(1)], fmiWindowFraction(1));
    elseif isGaussianSkew
        fmiWindowUniform1D = gaussianWindowFilter_asymmetric(numRings(1), fmiWindowFraction);
        fmiWindowUniform = repmat(fmiWindowUniform1D, numWedges(1), 1);
    else
        fmiWindowUniform = gaussianWindowFilter([numWedges(1), numRings(1)], fmiWindowFraction, 'fraction');
    end
    
    % Create the spectral energy filter for the FMI images
    fmiSpectralFilterUniform = spectralEnergyFilter(numWedges(1), numRings(1), fmiRPCDiameter); 
else
    fmiWindowUniform = 0;
    fmiSpectralFilterUniform = 0;  
end

% Initialize peak height ratio vector
fmcPeakRatio = zeros(nImages, 1);

% Initialize the RPC peak height ratio vector
rpcPeakRatio = zeros(nImages, 1);

% Set the fmc difference method to 2,
% which means "forward difference"
% since this is a Lagrangian
% Monte Carlo test.
fmcDifferenceMethod = 2;

% Do the processing
if NPROCESSORS > 1

    % Loop through images
    parfor k = 1: nImages
        
       % Read in the first image in the pair
        IMAGE1 = imageMatrix1(:, :, k);

        % Read in the second image in the pair
        IMAGE2 = imageMatrix2(:, :, k);

        % Make the windows if needed
        if uniformWindowDistribution == 1
            fmiWindow = fmiWindowUniform;
            fmiSpectralFilter = fmiSpectralFilterUniform;
        else
            fmiWindow = gaussianWindowFilter([numWedges(k), numRings(k)], fmiWindowFraction, 'fraction');
            fmiSpectralFilter = spectralEnergyFilter(numWedges(k), numRings(k), fmiRPCDiameter);
        end

        % Do the FMC correlation
        [estimatedTranslationY(k), estimatedTranslationX(k),...
        estimatedRotation(k), estimatedScaling(k), ...
        fmcPeakRatio(k), rpcPeakRatio(k), ~] = ...
        ...
        FMC(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
        fmiWindow, fmiSpectralFilter, ...
        xImage, yImage, ...
        spectrum_height, spectrum_width,...
        xLP, yLP, rMin, rMax, fmcDifferenceMethod, run_compiled);
        
    end % End (parfor k = 1 : nImages )    

else % Else if nProcessors == 1
   
    for k = 1 : nImages
        % Read in the first image in the pair
        IMAGE1 = imageMatrix1(:, :, k);

        % Read in the second image in the pair
        IMAGE2 = imageMatrix2(:, :, k);

        % Make the windows if needed
        if uniformWindowDistribution == 1
            fmiWindow = fmiWindowUniform;
            fmiSpectralFilter = fmiSpectralFilterUniform;
        else
            fmiWindow = gaussianWindowFilter([numWedges(k), numRings(k)], fmiWindowFraction, 'fraction');
            fmiSpectralFilter = spectralEnergyFilter(numWedges(k), numRings(k), fmiRPCDiameter);
        end
        
        % Do the FMC correlation
        [estimatedTranslationY(k), estimatedTranslationX(k),...
        estimatedRotation(k), estimatedScaling(k), ...
        fmcPeakRatio(k), rpcPeakRatio(k), ~] = ...
        ...
        FMC(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
        fmiWindow, fmiSpectralFilter, ...
        xImage, yImage, ...
        spectrum_height, spectrum_width,...
        xLP, yLP, rMin, rMax, fmcDifferenceMethod, run_compiled);   

    end % End (for k = 1 : nImages)

end % End ( if nProcessors > 1)

rotationError = angleAbsDiff(trueRotation, estimatedRotation);
scalingError = trueScaling - estimatedScaling;
translationErrorX = trueTranslationX - estimatedTranslationX;
translationErrorY = trueTranslationY - estimatedTranslationY;

fxTheory = 1 * (1 - numRings(1)) * log(trueScaling) / log(rMax / rMin);

particleDiameter = Parameters.ParticleDiameterMean;


% % UNCOMMENT THIS WHEN NOT RUNNING COMPILED CODES
% save(SAVEPATH, 'JobFile', 'concentration', 'trueRotation', 'shearX', 'shearY', ...
%     'trueScaling', 'trueTranslationX', 'trueTranslationY', ...
%   'estimatedRotation', ...
% 'estimatedScaling', 'estimatedTranslationX', 'estimatedTranslationY',  'rotationError', ...
%  'scalingError', 'translationErrorX', 'translationErrorY',  ...
%  'fmcPeakRatio', 'imageHeight', 'imageWidth', ...
% 'rMin', 'rMax', 'numWedges', 'numRings', 'particleDiameter', 'fxTheory');

save(SAVEPATH, 'JobFile', 'trueRotation', ...
    'trueScaling', 'trueTranslationX', 'trueTranslationY', ...
  'estimatedRotation', ...
'estimatedScaling', 'estimatedTranslationX', 'estimatedTranslationY',  'rotationError', ...
 'scalingError', 'translationErrorX', 'translationErrorY', 'imageHeight', 'imageWidth', ...
'rMin', 'rMax', 'numWedges', 'numRings', 'particleDiameter', 'fxTheory');



end




