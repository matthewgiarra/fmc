% function runFmiErrorStatistics(SETTYPE, DIMENSIONS, STARTSET, ENDSET)
function runFmiErrorStatistics(JOBLIST)

nJobs = length(JOBLIST);

for n = 1 : nJobs
    
JOBFILE = JOBLIST(n);
flipYTranslation = JOBFILE.JobOptions.FlipYTranslation;

% Determine image height and width (for naming purposes)
regionHeight = JOBFILE.Parameters.RegionHeight;
regionWidth = JOBFILE.Parameters.RegionWidth;
imageType = JOBFILE.ImageType;
setType = JOBFILE.SetType;
caseName = JOBFILE.CaseName;
correlationType = JOBFILE.CorrelationType;
startSet = JOBFILE.Parameters.Sets.Start;
endSet = JOBFILE.Parameters.Sets.End;

fprintf(1, 'Case %s\n', caseName);

% Number of digits in file names
setDigits = 5;

% Main code repository 
Repository = determineLocalRepositoryPath;

% String format of dataset titles
setFormat = ['%0' num2str(setDigits) '.0f'];

% Directory containing data
dataDirectory = fullfile(Repository, 'analysis', 'data', imageType, setType, caseName,  [num2str(regionHeight) 'x' num2str(regionWidth)], correlationType);

% Write directory
writeDir = fullfile(Repository, 'analysis', 'data',  imageType, setType, caseName, [num2str(regionHeight) 'x' num2str(regionWidth)], 'stat');
if ~exist(writeDir, 'dir')
    mkdir(writeDir)
end

% Number of data sets
nSets = endSet - startSet + 1; % Number of sets

% First file in the set
firstFile = fullfile( dataDirectory, ['errorAnalysis_' setType '_' correlationType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_' num2str(1, setFormat) '.mat'] ); % Name of first file
    
if exist(firstFile, 'file') % Load the file if it exists
    load(firstFile)
else
    fprintf(1, 'File does not exist!\n');
end
    
% Number of images is the number of elements in any of the results vectors
nImages = length(trueScaling); 

% % Initialize matrices to hold error analysis data 

% True rotations
trueRotationMat = zeros(nImages, nSets);
trueRotationPixelsMat = zeros(nImages, nSets);

% True scalings
trueScalingMat = zeros(nImages, nSets);
trueScalingPixelsMat = zeros(nImages, nSets);

% True translations
trueTranslationXMat = zeros(nImages, nSets);
trueTranslationYMat = zeros(nImages, nSets);

% true shear
trueShearXMat = zeros(nImages, nSets);
trueShearYMat = zeros(nImages, nSets);

% Estimated translations
estimatedTranslationXMat = zeros(nImages, nSets);
estimatedTranslationYMat = zeros(nImages, nSets);

% Estimated rotation
estimatedRotationMat = zeros(nImages, nSets);
estimatedScalingMat = zeros(nImages, nSets);

% Translation Error
translationAbsErrorXMat = zeros(nImages, nSets);
translationAbsErrorYMat = zeros(nImages, nSets);

% Rotation error
rotationAbsErrorMat = zeros(nImages, nSets);
rotationPixelsAbsErrorMat = zeros(nImages, nSets);

% Scaling error
scalingAbsErrorMat = zeros(nImages, nSets);
scalingPixelsAbsErrorMat = zeros(nImages, nSets);
 
 % Particle concentration
concentrationMat = zeros(nImages, nSets);

% Common signal
NiFiFsMat = zeros(nImages, nSets);

% Signal to noise ratios
peakRatioMat = zeros(nImages, nSets);

% Peak shifts in the FMI Correlation Plane
fmiTranslationXMat = zeros(nImages, nSets);
fmiTranslationYMat = zeros(nImages, nSets);
numWedgesMat = zeros(nImages, nSets);
numRingsMat = zeros(nImages, nSets);
spatialWindowFractionMat = zeros(2, nSets);
fmiWindowFractionMat = zeros(2, nSets);
spatialRPCdiameterMat = zeros(1, nSets);
fmiRPCdiameterMat = zeros(1, nSets);
imageWidth = zeros(1, nSets);
imageHeight = zeros(1, nSets);


for k = startSet : endSet
    
% Inform the user
    if rem(k, 100) == 0
        fprintf(1, ['Analyzing set ' setType '_' correlationType '_' num2str(k, setFormat) '\n' ], k);
    end
   
% Speify path to data
     dataPath = fullfile( dataDirectory, [ 'errorAnalysis_' setType '_' correlationType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_' num2str( k, setFormat ) '.mat' ] ); 
        
% Load the data
    Data = load(dataPath);
            
% True translations
    trueTranslationXMat( : , k ) = Data.trueTranslationX;
    
% True Y translation. Flip direction if specified.
    if flipYTranslation
        trueTranslationYMat( : , k ) = -1 * Data.trueTranslationY;
    else
        trueTranslationYMat( : , k ) =      Data.trueTranslationY;
    end
        
        
% True rotations
    trueRotationMat( : , k ) = Data.trueRotation;

% True scalings
    trueScalingMat( : , k ) = Data.trueScaling;
    
% True shearing
    trueShearXMat(:, k) = Data.shearX;
    trueShearYMat(:, k) = Data.shearY;
            
% Estimated translations
    estimatedTranslationXMat( : , k)  = Data.estimatedTranslationX;
    estimatedTranslationYMat( : , k ) = Data.estimatedTranslationY;
      
% Particle concentrations
    concentrationMat( : , k ) = Data.concentration;
        
%  Peak Ratio
%     peakRatioMat ( : , k) = Data.peakRatio;
    
% Processing parameters   
    imageHeight( k ) = Data.imageHeight;
    imageWidth( k ) = Data.imageWidth;
      
    if max(Data.translationErrorX) > imageHeight(k)
        keyboard
    end
    
% Average number of particle per windowed interrogation region
%     Ni =  (Data.concentration .* Data.imageHeight .* Data.imageWidth) * Data.spatialWindowFraction( 1 ) * Data.spatialWindowFraction( 2 );

% Loss of signal due to horizontal translation
%     Fix = 1 - abs(Data.trueTranslationX) ./ (Data.imageWidth * Data.spatialWindowFraction(2) ); % Loss of correlation due to horizontal in-plane motion

% Loss of signal due to vertical translation        
%     Fiy = 1 - abs( Data.trueTranslationY ) ./ (Data.imageHeight * Data.spatialWindowFraction(1) ); % Loss of correlation due to vertical in-plane motion

% Loss of signal due to scaling
%     Fs = 1 - abs(1 - Data.trueScaling) ./ (1 + Data.trueScaling); % Loss of correlation due to scaling
        
% Total estimated common signal
%     NiFiFsMat( :, k) = Ni .* Fix .* Fiy .* Fs;
   
% Deal with the FMC-specific parameters
    if strcmpi(correlationType, 'fmc');

        % Estimated rotations
        estimatedRotationMat( : , k ) = Data.estimatedRotation;
        
        % Estimated scaling
        estimatedScalingMat( : , k ) = Data.estimatedScaling;

        % Rotation errors
        rotationAbsErrorMat( : , k ) = abs(Data.rotationError);

        % Scaling errors
        scalingAbsErrorMat( : , k ) = abs(Data.scalingError);

        % Numbers of rings and wedges
        numRingsMat(:,  k ) = Data.numRings;
        numWedgesMat( :, k ) = Data.numWedges;

    end

end % End (For k = STARTSET : ENDSET )
    
% % Reshape the compiled statistics into vectors

fprintf(1, 'Done. \n \n');
% Translation errors.
translationAbsErrorXMat = abs(trueTranslationXMat - estimatedTranslationXMat);
translationAbsErrorYMat = abs(trueTranslationYMat - estimatedTranslationYMat);

% True horizontal translation 
trueTranslationX = reshape(trueTranslationXMat, numel(trueTranslationXMat), 1);

% True vertical translation
trueTranslationY = reshape(trueTranslationYMat, numel(trueTranslationYMat), 1);

% True Rotation
trueRotation = reshape(trueRotationMat, numel(trueRotationMat), 1);

% True rotation (pixels)
trueRotationPixels = reshape(trueRotationPixelsMat, numel(trueRotationPixelsMat), 1);

% True Scaling
trueScaling = reshape(trueScalingMat, numel(trueScalingMat), 1);

% True shearing
trueShearX = reshape(trueShearXMat, numel(trueShearXMat), 1);
trueShearY = reshape(trueShearYMat, numel(trueShearYMat), 1);

% Estimated horizontal translation
estimatedTranslationX = reshape(estimatedTranslationXMat, numel(estimatedTranslationXMat), 1);

% Estimated vertical translation
estimatedTranslationY = reshape(estimatedTranslationYMat, numel(estimatedTranslationYMat), 1);

% Error in horizontal translation estimate
translationAbsErrorX = reshape(translationAbsErrorXMat, numel(translationAbsErrorXMat), 1);

% Error in vertical translation estimate
translationAbsErrorY = reshape(translationAbsErrorYMat, numel(translationAbsErrorYMat), 1);

% Peak Ratio
% peakRatio = reshape(peakRatioMat, numel(peakRatioMat), 1);

% Common signal between images
% NiFiFs = reshape(NiFiFsMat, numel(NiFiFsMat), 1);

% Particle concentration
concentration = reshape(concentrationMat, numel(concentrationMat), 1);

% Mean true translation
trueTranslationXMean = (mean(trueTranslationXMat, 1))';
trueTranslationYMean = (mean(trueTranslationYMat, 1))';

% Mean true rotations
trueRotationMean = (mean(trueRotationMat, 1))';
% trueRotationPixelsMean = (mean(trueRotationPixelsMat, 1))';

% Mean true scalings
trueScalingMean = (mean(trueScalingMat, 1))';
% trueScalingPixelsMean = (mean(trueScalingPixelsMat, 1))';

% Deal with FMC-specific stuff
if strcmpi(correlationType, 'fmc')

    % Estimated Rotation
    estimatedRotation = reshape(estimatedRotationMat, numel(estimatedRotationMat), 1);

    % Estimated scaling
    estimatedScaling = reshape(estimatedScalingMat, numel(estimatedScalingMat), 1);

    % Error in rotation estimate
    rotationAbsError = reshape(rotationAbsErrorMat, numel(rotationAbsErrorMat), 1);

    % Error in rotation estimate (pixels)
    rotationPixelsAbsError = reshape(rotationPixelsAbsErrorMat, numel(rotationPixelsAbsErrorMat), 1);

    % Error in scaling estimate
    scalingAbsError = reshape(scalingAbsErrorMat, numel(scalingAbsErrorMat), 1);

    % Number of wedges
    numWedges = reshape(numWedgesMat, numel(numWedgesMat), 1);
    numRings =  reshape(numRingsMat, numel(numRingsMat), 1); 
    
    % I'm not sure what I had in mind when I calculated all of these mean and
    % std values below. I decided against saving these variables to disk for now but left
    % the code intact in case I decide to use them later.
 
%     % Mean and std-deviation FMI correlation peak shifts
%     fmiTranslationXMean = (mean(fmiTranslationXMat, 1))';
%     fmiTranslationXStd = (std(fmiTranslationXMat, 1))';
%     fmiTranslationYMean = (mean(fmiTranslationYMat, 1))';
%     fmiTranslationYStd = (std(fmiTranslationYMat, 1))';
% 
%     % Mean and std-deviation of estimated rotations
%     estimatedRotationMean = (mean(estimatedRotationMat, 1))';
%     estimatedRotationStd = (std(estimatedRotationMat, 1))';

%     % Mean and std-deviation of rotation transformation errors
%     rotationAbsErrorMean = (mean(rotationAbsErrorMat, 1))';
%     rotationAbsErrorStd = (std(rotationAbsErrorMat, 1))';
% 
%     % Mean and std-deviations of rotation pixel-displacement errors
%     rotationPixelsAbsErrorMean = (mean(rotationPixelsAbsErrorMat, 1))';
%     rotationPixelsAbsErrorStd = (std(rotationPixelsAbsErrorMat, 1))';
% 
%     % Mean and std-deviation of estimated scalings
%     estimatedScalingMean = (mean(estimatedScalingMat, 1))';
%     estimatedScalingStd = (std(estimatedScalingMat, 1))';
% 
%     % Mean and std-deviation of scaling errors
%     scalingAbsErrorMean = (mean(scalingAbsErrorMat, 1))';
%     scalingAbsErrorStd = (std(scalingAbsErrorMat, 1))';
% 
%     % Mean and std-deviations of scaling pixel-displacement errors
%     scalingPixelsAbsErrorMean = (mean(scalingPixelsAbsErrorMat, 1))';
%     scalingPixelsAbsErrorStd = (std(scalingPixelsAbsErrorMat, 1))';
%     
end

Ext = '.mat';
archiveDir = fullfile(writeDir, 'archive');
saveFileName = [ 'errorStatistics_' setType '_'...
    correlationType '_h' num2str(regionHeight) '_w' num2str(regionWidth)...
    '_' num2str(startSet, setFormat) '-' num2str(endSet, setFormat) ];
saveFilePath = fullfile(writeDir, [saveFileName Ext]);

% Archive any existing files
if exist(saveFilePath, 'file')
    if ~exist(archiveDir, 'dir')
        mkdir(archiveDir)
    end;
    
   dateString = datestr(now);
   dateString = strrep(dateString, ' ', '_');
   dateString = strrep(dateString, ':', '-');
   archiveFileName = [saveFileName '_' dateString Ext ];
   archiveFilePath = fullfile(writeDir, 'archive', archiveFileName);

   movefile(saveFilePath, archiveFilePath);
    
end

if strcmpi(correlationType, 'fmc')

 save(saveFilePath,...
    'trueTranslationX', 'trueTranslationY', ...
    'trueRotation', 'trueScaling', ...
    'estimatedTranslationX', 'estimatedTranslationY', ...
    'estimatedScaling', 'estimatedRotation', ...
    'translationAbsErrorX', 'translationAbsErrorY', ...
    'rotationAbsError', 'scalingAbsError', ...
     'imageHeight', 'imageWidth', 'concentration', 'peakRatio', 'numWedges', ...
     'numRings', 'trueShearX', 'trueShearY');
 
else % Save error statistics (RPC or other non-FMC)
    save(saveFilePath,...
    'trueTranslationX', 'trueTranslationY', ...
    'trueRotation', 'trueScaling', ...
    'estimatedTranslationX', 'estimatedTranslationY', ...
    'translationAbsErrorX', 'translationAbsErrorY', ...
     'imageHeight', 'imageWidth', 'concentration', ...
     'trueRotationPixels', 'trueShearX', 'trueShearY');
 
end

end

end % END OF FUNCTION










