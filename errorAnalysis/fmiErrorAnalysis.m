function fmiErrorAnalysis(IMDIR, SAVEPATH, IMBASE, NDIGITS, EXTENSION, PARAMETERSPATH, SPATIALWINDOW, FMIWINDOW, IMAGESPECTRALFILTER, FMISPECTRALFILTER, CONSTANTS, SAVEDATA, NPROCESSORS)

% Default to not saving data
if nargin < 11
    SAVEDATA = 0;
end

% Input parameters that don't change
% SPATIALRPCDIAMETER = 2.8; % Spatial image RPC diameter (pixels)
% SPATIALWINDOWFRACTION = [0.5 0.5]; % Spatial image window fraction (y, x)
% FMIRPCDIAMETER = 3.3; % FMI RPC Diameter (pixels)
% FMIWINDOWFRACTION = [1 0.5]; % FMI Window Fraction (y, x)
% LSSWINDOWSCALE = 20; % Tuckey-window scaling for the Least Squares Solver images

% Specify numbering format based on number of digits
numberFormat = ['%0' num2str(NDIGITS) '.0f'];

% try
    
    % Load in known transformation parameters (used for error analysis)
    Parameters = load(PARAMETERSPATH);
    trueRotation = Parameters.Rotation;
    trueScaling = Parameters.Scaling;
    concentration = Parameters.Concentration;
    trueTransformationMatrix = Parameters.Tforms;
    trueTranslationX = Parameters.TranslationX;
    trueTranslationY = Parameters.TranslationY;
    % Size of raw images
    imageHeight = Parameters.ImageHeight;
    imageWidth = Parameters.ImageWidth;

    % Number of images
    nImages = length(trueRotation); 

    % Start and stop files
    startFile = 1;
    stopFile = nImages;

    % Paths to input images (first images in the pairs)
    REFERENCEFILEPATH = fullfile(IMDIR, 'referenceImage.tiff');  
    
    % Paths to input images (second images in the pairs)
    FILEPATHSB = makeFilepaths(IMDIR, IMBASE, numberFormat, EXTENSION, startFile, stopFile, '', 1);

    % Radius of largest circle that can be inscribed on the images. This is
    % used for the FMI transformation.
    rMax = min(imageHeight / 2, imageWidth / 2);

    % Initialize vectors to hold estimates
    estimatedRotationInitial = zeros(nImages, 1);
    estimatedRotation = zeros(nImages, 1);
    estimatedScaling = zeros(nImages, 1);
    estimatedTranslationX = zeros(nImages, 1);
    estimatedTranslationY = zeros(nImages, 1);

    % Initialize vectors to hold raw translation estimates
    fmiTranslationY = zeros(nImages, 1);
    fmiTranslationX = zeros(nImages, 1);

    % Initialize vectors to hold true translations
    trueRotationPixels = zeros(nImages, 1);
    trueScalingPixels = zeros(nImages, 1);

    % Initialize vectors to hold absolute transformation errors
    rotationAbsError = zeros(nImages, 1);
    scalingAbsError = zeros(nImages, 1);
    translationAbsErrorX = zeros(nImages, 1);
    translationAbsErrorY = zeros(nImages, 1);

    % Initialize vectors to hold absolute pixel errors
    rotationPixelsAbsError = zeros(nImages, 1);
    scalingPixelsAbsError = zeros(nImages, 1);
    
    % Initialize vector to hold peak ratios
    peakRatio = zeros(nImages, 1);
    
%     % Initialize transformation structure
%     TFORM = maketform('affine',  [ 1 0 0 ; 0 1 0 ; 0 0 1 ] );
    
    if NPROCESSORS > 1

        matlabpool( NPROCESSORS ) % Open a matlab pool

        % Read in the first image in the pair
        IMAGE1 = imread( REFERENCEFILEPATH );    

        % Convert first image to FMI   
        FMI1 = FMIWINDOW .* im2fmi( SPATIALWINDOW .* double( IMAGE1 ) );
        
        % Size of FMI images
        [fmiHeight fmiWidth] = size(FMI);
        
        % Loop through images
        parfor k = 1: nImages

            if rem( k, 10 ) == 0
                fprintf( 1, 'Analyzing image %6.0f of %6.0f\n', k, nImages );
            end

        % Read in the second image in the pair
            IMAGE2 = imread(FILEPATHSB(k, :)); 

        % Convert second image to FMI
            FMI2 = FMIWINDOW .* im2fmi(SPATIALWINDOW .* double(IMAGE2));

        % Calculate FMI RPC between images     
            [estimatedRotationInitial(k) estimatedScaling(k) fmiTranslationY(k) fmiTranslationX(k) correlationPlane] = FMIRPC(FMI1, FMI2, FMISPECTRALFILTER);

        % Resolve rotation ambiguity
             [estimatedRotation(k) estimatedTranslationX(k) estimatedTranslationY(k) peakRatio(k)] = resolveRotationAmbiguity(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER, estimatedRotationInitial(k), estimatedScaling(k));

        % Calculate absolute errors (i.e. absolute difference between the calculated and true values for scaling and rotation)
            rotationAbsError(k) = abs(trueRotation(k) - estimatedRotation(k));
            scalingAbsError(k) = abs(trueScaling(k) - estimatedScaling(k));
            translationAbsErrorX(k) = abs( trueTranslationX(k)  - estimatedTranslationX(k)) ;
            translationAbsErrorY(k) = abs( trueTranslationY(k) - estimatedTranslationY(k)) ;

        % Calculated absolute errors of pixel shift in FMI correlation plane
        % between true and calculated values of scaling and rotation
            trueRotationPixels(k) = trueRotation(k) * fmiHeight / (2 * pi); % Theoretical peak displacement due to true rotation
            rotationPixelsAbsError(k) = abs(trueRotationPixels(k) - fmiTranslationY(k)); % Absolute difference between theoretical and measured peak displacements
            trueScalingPixels(k) = - fmiWidth * log(trueScaling(k)) / log(rMax); % Theoretical peak displacement due to true scaling
            scalingPixelsAbsError(k) = abs(trueScalingPixels(k) - fmiTranslationX(k)); % Absolute difference between theoretical and measured peak displacements

        end % End (parfor k = 1 : nImages )
        
        matlabpool close;

        else % Else if nProcessors == 1

        % Read in the first image in the pair
        IMAGE1 = imread(REFERENCEFILEPATH);  
           
        % Convert first image to FMI   
        FMI1 = FMIWINDOW .* im2fmi(SPATIALWINDOW .* double(IMAGE1));
        
       % Size of FMI images
        [fmiHeight fmiWidth] = size(FMI1);
        
        for k = 1: nImages
                if rem(k, 10) == 0  
                    fprintf(1, 'Analyzing image %6.0f of %6.0f\n', k, nImages);
                end

    % Read in the second image in the pair
                IMAGE2 = imread(FILEPATHSB(k, :)); 

    % Convert second image to FMI
                FMI2 = FMIWINDOW .* im2fmi(SPATIALWINDOW .* double(IMAGE2));

    % Calculate FMI RPC between images     
              [estimatedRotationInitial(k) estimatedScaling(k) fmiTranslationY(k) fmiTranslationX(k) correlationPlane] = FMIRPC(FMI1, FMI2, FMISPECTRALFILTER);

    % Resolve rotation ambiguity
              [estimatedRotation(k) estimatedTranslationX(k) estimatedTranslationY(k) peakRatio(k) IMAGEOUT] = resolveRotationAmbiguity(IMAGE1, IMAGE2, SPATIALWINDOW, IMAGESPECTRALFILTER, estimatedRotationInitial(k), estimatedScaling(k));

%               
%             compositeImage = uint8(cat(3, IMAGE2, IMAGEOUT, zeros(size(IMAGE2))));
%             imshow(compositeImage);
%             set(gcf, 'OuterPosition', [1 1 800 800]);
%             print(1, '-opengl', '-dpng', '-r500',  ['~/Desktop/images/image_' num2str(k, '%04.0f') '.png']);
             
              
    % Calculate absolute errors (i.e. absolute difference between the calculated and true values for scaling and rotation)
                rotationAbsError(k) = abs(trueRotation(k) - estimatedRotation(k));
                scalingAbsError(k) = abs(trueScaling(k) - estimatedScaling(k));
                translationAbsErrorX(k) = abs( trueTranslationX(k)  - estimatedTranslationX(k)) ;
                translationAbsErrorY(k) = abs( trueTranslationY(k) - estimatedTranslationY(k)) ;

    % Calculated absolute errors of pixel shift in FMI correlation plane
    % between true and calculated values of scaling and rotation
                trueRotationPixels(k) = trueRotation(k) * fmiHeight / (2 * pi); % Theoretical peak displacement due to true rotation
                rotationPixelsAbsError(k) = abs(trueRotationPixels(k) - fmiTranslationY(k)); % Absolute difference between theoretical and measured peak displacements
                trueScalingPixels(k) = - fmiWidth * log(trueScaling(k)) / log(rMax); % Theoretical peak displacement due to true scaling
                scalingPixelsAbsError(k) = abs(trueScalingPixels(k) - fmiTranslationX(k)); % Absolute difference between theoretical and measured peak displacements

        end % End (for k = 1 : nImages)

    end % End ( if nProcessors > 1)
    
    % Extract constants from structure
    spatialWindowFraction = CONSTANTS.SpatialWindowFraction;
    fmiWindowFraction = CONSTANTS.FmiWindowFraction;
    spatialRPCdiameter = CONSTANTS.SpatialRPCdiameter;
    fmiRPCdiameter = CONSTANTS.FmiRPCdiameter;

    % rename image directory variable
    imageDirectory = IMDIR;
    
if SAVEDATA
    % Save outputs to disk
    save( SAVEPATH, 'imageDirectory', 'concentration', 'trueRotation', 'trueRotationPixels', 'trueScaling', 'trueTranslationX', 'trueTranslationY', 'trueScalingPixels', 'fmiTranslationY', 'fmiTranslationX', 'estimatedRotation', 'estimatedScaling', 'estimatedTranslationX', 'estimatedTranslationY',  'rotationAbsError', 'rotationPixelsAbsError', 'scalingAbsError', 'scalingPixelsAbsError', 'translationAbsErrorX', 'translationAbsErrorY',  'peakRatio', 'spatialWindowFraction', 'fmiWindowFraction', 'spatialRPCdiameter', 'fmiRPCdiameter', 'imageHeight', 'imageWidth', 'fmiHeight', 'fmiWidth');
end

% catch err;
%     
% % Open error log file and write error message
% fclose all;
% fid = fopen(fullfile( '~/Desktop', 'errorLog.txt' ), 'a' );
% fprintf( fid, ['Problem with set ' IMDIR '\n' err.message '\n\n'] );
% fclose(fid);
% end



end




