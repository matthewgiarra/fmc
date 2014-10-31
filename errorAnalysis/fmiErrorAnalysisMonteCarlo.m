function fmiErrorAnalysisMonteCarlo(JOBFILE, SAVEPATH, IMAGEFILEPATH, NDIGITS, PARAMETERSPATH, SPATIALWINDOW, IMAGESPECTRALFILTER, CONSTANTS, NPROCESSORS)

% Rename the Jobfile variable.
JobFile = JOBFILE;

% Flag specifying whether or not to run the nonlinear least square
% solver and find the full affine transformation.
doAffineTransform = JobFile.JobOptions.DoAffineTransform;

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

% Noise mean and standard deviation fractions
noiseMean = CONSTANTS.NoiseMean;
noiseStd = CONSTANTS.NoiseStd;

% Specify numbering format based on number of digits
numberFormat = ['%0' num2str(NDIGITS) '.0f'];

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

% Number of images
nImages = length(trueRotation); 

% Start and stop files
startFile = 1;
stopFile = nImages;

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
resampledSourceHeight = fftSize(1);
resampledSourceWidth  = fftSize(2);


% Minimum and maximum radii of log-polar resampling of the FFT magnitude. 
rMax = min(resampledSourceHeight, resampledSourceWidth) / 2 - 1;
rMin = CONSTANTS.RMin;

% Make the image coordinates
[xImage, yImage] = meshgrid(1:imageWidth, 1:imageHeight);

% FFT Spectrum coordinates
[resampledSourceX, resampledSourceY] = meshgrid(1:resampledSourceWidth, 1:resampledSourceHeight);
% [rsX, rsY] = meshgrid(1:resampledSourceWidth, 1:resampledSourceHeight);
% 
% resampledSourceX = fftshift(rsX);
% resampledSourceY = fftshift(rsY);

% Figure out the log polar resampling coordinates.
[xLP, yLP] = LogPolarCoordinates([resampledSourceHeight, resampledSourceWidth], numWedges(1), numRings(1), rMin, rMax, 2*pi);
% [xlp, ylp] = LogPolarCoordinates([resampledSourceHeight, resampledSourceWidth], numWedges(1), numRings(1), rMin, rMax, pi);
% xLP = fftshift(xlp, 2);
% yLP = (ylp);

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
isGaussianSym = ~isempty(regexpi(fmiWindowType, 'gauss_sym'));

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
    
%     fmiWindowUniform = gaussianWindowFilter([numWedges(1), numRings(1)], fmiWindowFraction, 'fraction');
    fmiSpectralFilterUniform = spectralEnergyFilter(numWedges(1), numRings(1), fmiRPCDiameter); 
else
    fmiWindowUniform = 0;
    fmiSpectralFilterUniform = 0;  
end

% Load the image matrix.
% imageFilePath = fullfile(IMDIR, ['raw_image_matrix_mc_h' num2str(imageHeight) '_w' num2str(imageWidth) '_seg_' num2str(1, numberFormat) '_' num2str(nImages, numberFormat) '.mat'] );
imageData = load(IMAGEFILEPATH);

% Matrices containing the image data for the first and second image in each pair.
imageMatrix1 = double(imageData.imageMatrix1);
imageMatrix2 = double(imageData.imageMatrix2);

% Max value of the images
maxVal = double(intmax(class(imageData.imageMatrix1)));

% Make noise matrices. The 2.8 corresponds to the multiple of the standard
% deviation corresponding to a 99.5% coverage factor.
% noiseMatrix1 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(imageMatrix1));
% noiseMatrix2 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(imageMatrix2));

% ftFilter = fftshift(makeFFTFilter([imageHeight, imageWidth], 1));
ftFilter = ones([numWedges(1), numRings(1)]);

% Initialize peak height ratio vector
fmcPeakRatio = zeros(nImages, 1);

% Initialize the RPC peak height ratio vector
rpcPeakRatio = zeros(nImages, 1);

% Initialize affine parameters matrix
affineParameters = zeros(nImages, 6);

RUN_COMPILED = ~ismac;

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

%         % Do the FMC correlation
%         [estimatedTranslationY(k), estimatedTranslationX(k),...
%         estimatedRotation(k), estimatedScaling(k), ...
%         fmcPeakRatio(k), rpcPeakRatio(k), ~] = ...
%         ...
%         FMC(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
%         fmiWindow, fmiSpectralFilter, ...
%         ftFilter, xImage, yImage, ...
%         resampledSourceX, resampledSourceY,...
%         xLP, yLP, rMin, rMax, doAffineTransform );
    
    
    
            % Do the FMC correlation
        [estimatedTranslationY(k), estimatedTranslationX(k),...
        estimatedRotation(k), estimatedScaling(k), ...
        fmcPeakRatio(k), rpcPeakRatio(k), ~] = ...
        ...
        FMC(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
        fmiWindow, fmiSpectralFilter, ...
        xImage, yImage, ...
        resampledSourceX, resampledSourceY,...
        xLP, yLP, rMin, rMax );
    
%     
%         [estimatedTranslationY(k), estimatedTranslationX(k),...
%         estimatedRotation(k), estimatedScaling(k)] = ...
%         FMC_compiled(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
%         fmiWindow, fmiSpectralFilter, ...
%         xImage, yImage, resampledSourceHeight, resampledSourceWidth,...
%         xLP, yLP, rMin, rMax, RUN_COMPILED ); 
%     
    
        
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
%         [estimatedTranslationY(k), estimatedTranslationX(k),...
%         estimatedRotation(k), estimatedScaling(k), ...
%         fmcPeakRatio(k), rpcPeakRatio(k), fmiPlane] = ...
%         FMC(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
%         fmiWindow, fmiSpectralFilter, ...
%         ftFilter, xImage, yImage, ...
%         resampledSourceX, resampledSourceY,...
%         xLP, yLP, rMin, rMax, doAffineTransform );   
   
            [estimatedTranslationY(k), estimatedTranslationX(k),...
        estimatedRotation(k), estimatedScaling(k), ...
        fmcPeakRatio(k), rpcPeakRatio(k)] = ...
        ...
        FMC(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
        fmiWindow, fmiSpectralFilter, ...
        xImage, yImage, ...
        resampledSourceX, resampledSourceY,...
        xLP, yLP, rMin, rMax, 2, 0);   
    
    
%     %       Do the FMC correlation
%         [estimatedTranslationY(k), estimatedTranslationX(k),...
%         estimatedRotation(k), estimatedScaling(k)] = ...
%         FMC_compiled(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER,...
%         fmiWindow, fmiSpectralFilter, ...
%         xImage, yImage, resampledSourceHeight, resampledSourceWidth,...
%         xLP, yLP, rMin, rMax, RUN_COMPILED);   
    
%         mesh(fmiPlane ./ max(fmiPlane(:)));
       
    end % End (for k = 1 : nImages)

end % End ( if nProcessors > 1)

rotationError = angleAbsDiff(trueRotation, estimatedRotation);
scalingError = trueScaling - estimatedScaling;
translationErrorX = trueTranslationX - estimatedTranslationX;
translationErrorY = trueTranslationY - estimatedTranslationY;

% fxTheory = 1 * (1 - numRings(1)) * log(trueScaling) / log(rMax / rMin);
fxTheory = 1 * (1 - numRings(1)) * log(trueScaling) / log(rMax / rMin);

particleDiameter = Parameters.ParticleDiameter(1);


% % UNCOMMENT THIS WHEN NOT RUNNING COMPILED CODES
save(SAVEPATH, 'JobFile', 'concentration', 'trueRotation', 'shearX', 'shearY', ...
    'trueScaling', 'trueTranslationX', 'trueTranslationY', ...
  'estimatedRotation', ...
'estimatedScaling', 'estimatedTranslationX', 'estimatedTranslationY',  'rotationError', ...
 'scalingError', 'translationErrorX', 'translationErrorY',  ...
 'fmcPeakRatio', 'imageHeight', 'imageWidth', ...
'rMin', 'rMax', 'numWedges', 'numRings', 'particleDiameter', 'fxTheory', 'affineParameters');

% save(SAVEPATH, 'JobFile', 'concentration', 'trueRotation', 'shearX', 'shearY', ...
%     'trueScaling', 'trueTranslationX', 'trueTranslationY', ...
%   'estimatedRotation', ...
% 'estimatedScaling', 'estimatedTranslationX', 'estimatedTranslationY',  'rotationError', ...
%  'scalingError', 'translationErrorX', 'translationErrorY', 'imageHeight', 'imageWidth', ...
% 'rMin', 'rMax', 'numWedges', 'numRings', 'particleDiameter', 'fxTheory', 'affineParameters');



end




