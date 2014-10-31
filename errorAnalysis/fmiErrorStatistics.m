setType = 'scaling';
height = 128; width = 128; % Image dimensions
conc = 0.050; % Particle concentration

nProcessors = 7;
startSet = 1;
endSet = 10000;
setDigits = 5;

% Main code repository
repository = fullfile(pwd, '..', '..'); 

setBase = [ setType '_' ];
IMBASE = [ setType '_' ];

setFormat = ['%0' num2str(setDigits) '.0f'];

dataDirectory = fullfile(repository, 'results', [num2str(height) 'x' num2str(width)], ['c_' num2str(conc, '%3.3f')], setType);

nSets = endSet - startSet + 1; % Number of sets

% Load in a data set and check for information
Found = 0; % Has the file been found?
Set = 1; % Set number

while Found == 0
    firstFile = fullfile( dataDirectory, ['errorAnalysis_' setType '_' num2str(height) 'x' num2str(width) '_c_' num2str(conc, '%3.3f') '_' num2str(Set, setFormat) '.mat'] ); % Name of first file
    
    if exist(firstFile, 'file') % Load the file if it exists
        load(firstFile)
        Found = 1;
    end
       
end
    
% Number of images is the number of elements in any of the results vectors
nImages = length(trueScaling); 

% Initialize matrices to hold error analysis data 
 trueRotationMat = zeros(nSets, nImages);
 trueRotationPixelsMat = zeros(nSets, nImages);
 
 trueScalingMat = zeros(nSets, nImages);
 trueScalingPixelsMat = zeros(nSets, nImages);
 
 trueShearXMat = zeros(nSets, nImages);
 trueShearYMat = zeros(nSets, nImages);
 
 fmiTranslationXMat = zeros(nSets, nImages);
 fmiTranslationYMat = zeros(nSets, nImages);
 
 estimatedRotationMat = zeros(nSets, nImages);
 estimatedScalingMat = zeros(nSets, nImages);
 
 rotationAbsErrorMat = zeros(nSets, nImages);
 rotationPixelsAbsErrorMat = zeros(nSets, nImages);
 
 scalingAbsErrorMat = zeros(nSets, nImages);
 scalingPixelsAbsErrorMat = zeros(nSets, nImages);
 
 spatialWindowFractionMat = zeros(nSets, 2);
 fmiWindowFractionMat = zeros(nSets, 2);
 spatialRPCdiameterMat = zeros(nSets, 1);
 fmiRPCdiameterMat = zeros(nSets, 1);
 imageWidth = zeros(nSets, 1);
 imageHeight = zeros(nSets, 1);
 
 
 
 matlabpool(nProcessors)
 
parfor k = startSet : endSet
    
    if rem(k, 100) == 0
        fprintf(1, ['Analyzing set ' setBase num2str(k, setFormat) '\n' ], k);
    end
   
    dataPath = fullfile( dataDirectory, [ 'errorAnalysis_' setBase num2str(height) 'x' num2str(width) '_c_' num2str(conc, '%3.3f') '_' num2str( k, setFormat ) '.mat' ] ); % Data path
    
    Data = load(dataPath);
    
    % True rotations
    trueRotationMat( k , : ) = Data.trueRotation;
    trueRotationPixelsMat( k, : ) = Data.trueRotationPixels;
    
    % True scalings
    trueScalingMat( k , : ) = Data.trueScaling;
    trueScalingPixelsMat( k, : ) = Data.trueScalingPixels;
    
    trueShearXMat(k, :) = Data.shearX;
    trueShearYMat(k, :) = Data.shearY;
 
    % FMI Correlation Peak Shifts
    fmiTranslationXMat( k, : ) = Data.fmiTranslationX;
    fmiTranslationYMat( k, : ) = Data.fmiTranslationY;
    
    % Estimated rotations and scalings
    estimatedRotationMat( k , : ) = Data.estimatedRotation;
    estimatedScalingMat( k , : ) = Data.estimatedScaling;
    
    % Rotation errors
    rotationAbsErrorMat( k , : ) = Data.rotationAbsError;
    rotationPixelsAbsErrorMat( k, : ) = Data.rotationPixelsAbsError;
    
    % Scaling errors
    scalingAbsErrorMat( k , : ) = Data.scalingAbsError;
    scalingPixelsAbsErrorMat( k, : ) = Data.scalingPixelsAbsError;
    
    spatialWindowFractionMat( k, : ) = Data.spatialWindowFraction;
    fmiWindowFractionMat( k, :) = Data.fmiWindowFraction;
    spatialRPCdiameterMat( k ) = Data.spatialRPCdiameter;
    fmiRPCdiameterMat( k ) = Data.fmiRPCdiameter;
    imageHeight( k ) = Data.imageHeight;
    imageWidth( k ) = Data.imageWidth;

end

matlabpool close

% Mean true rotations
trueRotationMean = (mean(trueRotationMat, 1))';
trueRotationPixelsMean = (mean(trueRotationPixelsMat, 1))';

% Mean true scalings
trueScalingMean = (mean(trueScalingMat, 1))';
trueScalingPixelsMean = (mean(trueScalingPixelsMat, 1))';

% Mean and std-deviation FMI correlation peak shifts
fmiTranslationXMean = (mean(fmiTranslationXMat, 1))';
fmiTranslationXStd = (std(fmiTranslationXMat, 1))';
fmiTranslationYMean = (mean(fmiTranslationYMat, 1))';
fmiTranslationYStd = (std(fmiTranslationYMat, 1))';

% Mean and std-deviation of estimated rotations
estimatedRotationMean = (mean(estimatedRotationMat, 1))';
estimatedRotationStd = (std(estimatedRotationMat, 1))';

% Mean and std-deviation of rotation transformation errors
rotationAbsErrorMean = (mean(rotationAbsErrorMat, 1))';
rotationAbsErrorStd = (std(rotationAbsErrorMat, 1))';

% Mean and std-deviations of rotation pixel-displacement errors
rotationPixelsAbsErrorMean = (mean(rotationPixelsAbsErrorMat, 1))';
rotationPixelsAbsErrorStd = (std(rotationPixelsAbsErrorMat, 1))';

% Mean and std-deviation of estimated scalings
estimatedScalingMean = (mean(estimatedScalingMat, 1))';
estimatedScalingStd = (std(estimatedScalingMat, 1))';

% Mean and std-deviation of scaling errors
scalingAbsErrorMean = (mean(scalingAbsErrorMat, 1))';
scalingAbsErrorStd = (std(scalingAbsErrorMat, 1))';

% Mean and std-deviations of scaling pixel-displacement errors
scalingPixelsAbsErrorMean = (mean(scalingPixelsAbsErrorMat, 1))';
scalingPixelsAbsErrorStd = (std(scalingPixelsAbsErrorMat, 1))';

% Save error statistics
save(fullfile(dataDirectory, [ 'errorStatistics_' setBase num2str(height) 'x' num2str(width) '_c_' num2str(conc, '%3.3f')  '.mat' ]),...
    'trueRotationMean', 'trueRotationPixelsMean', ...
    'trueScalingMean', 'trueScalingPixelsMean', ...
    'fmiTranslationXMean', 'fmiTranslationYMean', ...
    'fmiTranslationXStd', 'fmiTranslationYStd', ...
    'estimatedRotationMean', 'estimatedRotationStd', ...  
    'rotationAbsErrorMean', 'rotationAbsErrorStd', ...
    'rotationPixelsAbsErrorMean', 'rotationPixelsAbsErrorStd', ...
     'estimatedScalingMean', 'estimatedScalingStd', ...
     'scalingAbsErrorMean', 'scalingAbsErrorStd', ...
     'scalingPixelsAbsErrorMean', 'scalingPixelsAbsErrorStd', ...
     'imageHeight', 'imageWidth');

% Save compiled error analyses
save(fullfile(dataDirectory, [ 'errorDataCompiled_' setBase num2str(height) 'x' num2str(width) '_c_' num2str(conc, '%3.3f')  '.mat' ]),...
  'trueRotationMat', 'trueRotationPixelsMat', ...
  'trueScalingMat', 'trueScalingPixelsMat', ...
  'fmiTranslationXMat', 'fmiTranslationYMat', ...
  'estimatedRotationMat', 'estimatedScalingMat', ...
  'rotationAbsErrorMat', 'scalingAbsErrorMat', ...
  'rotationPixelsAbsErrorMat', 'scalingPixelsAbsErrorMat', ...
  'spatialWindowFractionMat', 'fmiWindowFractionMat', ...
  'spatialRPCdiameterMat', 'fmiRPCdiameterMat', ...
  'imageHeight', 'imageWidth');










