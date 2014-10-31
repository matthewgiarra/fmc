% function runFmiJobFileMc(SETTYPE, CASENAME, DIMENSIONS, STARTSET, ENDSET, NPROCESSORS)
function runFmcJobFileMc(JOBLIST)

nJobs = length(JOBLIST);

for n = 1 : nJobs
    
JOBFILE = JOBLIST(n);

SAVEDATA = 1; % Force data save (hard code so that I don't spend a week processing data and forget to save the results)

skipExistingSets = JOBFILE.JobOptions.SkipExistingSets;
regionHeight = JOBFILE.Parameters.RegionHeight;
regionWidth = JOBFILE.Parameters.RegionWidth;
caseName = JOBFILE.CaseName;
correlationType = JOBFILE.CorrelationType;
imageType = JOBFILE.ImageType;
setType = JOBFILE.SetType;
startSet = JOBFILE.Parameters.Sets.Start;
endSet = JOBFILE.Parameters.Sets.End;
imagesPerSet = JOBFILE.Parameters.Sets.ImagesPerSet;
nProcessors = JOBFILE.JobOptions.NumberOfProcessors;

if JOBFILE.JobOptions.RepositoryPathIsAbsolute
    Repository = JOBFILE.Parameters.RepositoryPath;
else
    % Main code repository
    Repository = determineLocalRepositoryPath;
end

% Case directory
caseDir = fullfile(Repository, 'analysis', 'data', imageType,setType, caseName);

% Image parent directory
imageParentDirectory = fullfile(caseDir, [num2str(regionHeight) 'x' num2str(regionWidth)], 'raw');

% Write directory
writeDir = fullfile(caseDir, [num2str(regionHeight) 'x' num2str(regionWidth)], correlationType);

% Make the write directory if it doesn't exist
if ~exist(writeDir, 'dir')
    mkdir(writeDir)
end

% Base names of image sets
setBase = [setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_'];

% Base names of results files
saveBase = ['errorAnalysis_' setType '_' correlationType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_'];

% Base names of images
IMBASE = ['img_' setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_'];

% Number of digits in the set names
setDigits = 5;
% Numbering format tag for image sets
setFormat = ['%0' num2str(setDigits) '.0f'];

% Number of digits in image names
imageDigits = 6;
imageNumberFormat = ['%0' num2str(imageDigits) '.0f'];
% IMAGESTEP = 1; % Not sure why this is specified. Image step should always be 1 for this analysis, I think.
EXTENSION = '.tiff';

spatialWindowType =  JOBFILE.Parameters.Processing.SpatialWindowType; % Spatial window type
spatialWindowFraction = JOBFILE.Parameters.Processing.SpatialWindowFraction; % Spatial image window fraction (y, x)

fmiWindowType = JOBFILE.Parameters.Processing.FMIWindowType; % FMI Window Type
fmiWindowFraction = JOBFILE.Parameters.Processing.FMIWindowFraction; % FMI Window Fraction (y, x)

spatialRPCdiameter = JOBFILE.Parameters.Processing.SpatialRPCDiameter; % Spatial image RPC diameter (pixels)
fmiRPCdiameter = JOBFILE.Parameters.Processing.FmiRPCDiameter; % FMI RPC Diameter (pixels)

% Make the spatial window
SPATIALWINDOW = gaussianWindowFilter( [regionHeight regionWidth], spatialWindowFraction, spatialWindowType );

% Make the image spectral filter
IMAGESPECTRALFILTER = spectralEnergyFilter(regionHeight, regionWidth, spatialRPCdiameter); % Raw image RPC spectral energy filter

% Write some constants to a structure
consts.SpatialWindowFraction = spatialWindowFraction;
consts.FmiWindowFraction = fmiWindowFraction;
consts.SpatialRPCdiameter = spatialRPCdiameter;
consts.FmiRPCdiameter = fmiRPCdiameter;
consts.NWmax =  JOBFILE.Parameters.Processing.NumberOfWedges.Max;
consts.NWmin =  JOBFILE.Parameters.Processing.NumberOfWedges.Min;
consts.NRmin =  JOBFILE.Parameters.Processing.NumberOfRings.Min;
consts.NRmax =  JOBFILE.Parameters.Processing.NumberOfRings.Max;
consts.FFTSize = JOBFILE.Parameters.Processing.FFTSize;
consts.RMin = JOBFILE.Parameters.Processing.MinimumRadius;
consts.FMIWindowType = JOBFILE.Parameters.Processing.FMIWindowType;
consts.NoiseMean = JOBFILE.Parameters.Processing.Noise.Mean;
consts.NoiseStd = JOBFILE.Parameters.Processing.Noise.Std;

% Close any open files
fclose all;

% Create error log file
fid = fopen(fullfile( '~/Desktop', 'FMIerrorLog.txt' ), 'w');

% Print header to error log file
fprintf( fid, [ 'FMI-RPC Processing Error Log \nCase: ' setBase '\n\n' ] ); 
fclose(fid);

% List of set numbers
setList = startSet : endSet;
nSets = length(setList);

for k = 1 : nSets             
    fprintf(1, ['Analyzing set ' ...
         correlationType ' ' caseName ' ' setBase num2str(setList(k), setFormat) ' (' num2str(k)...
         ' of ' num2str(nSets) ')... ']); % Display status

     IMDIR = fullfile( imageParentDirectory, [ setBase num2str( setList(k), setFormat ) ], 'raw' ); % Image directory
     SAVEPATH = fullfile( writeDir, [ saveBase num2str( setList(k), setFormat ) '.mat' ] ); % Save path

    % Path to image parameters
    parametersPath = fullfile(imageParentDirectory,...
            [setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_' num2str(setList(k), setFormat)], 'parameters', ...
            ['imageParameters_' setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_seg_' num2str(1, imageNumberFormat) '_' num2str(imagesPerSet, imageNumberFormat) '.mat']);                         
    
    % Path to the raw image file
    imageFilePath = fullfile(IMDIR, ['raw_image_matrix_' setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_seg_' num2str(1, imageNumberFormat) '_' num2str(imagesPerSet, imageNumberFormat) '.mat'] );
        
        % Perform analysis on the image set if "skip existing sets" isn't
        % selected or if the set results don't exist.
    if (~skipExistingSets) || (skipExistingSets && ~exist(SAVEPATH, 'file'))

        % Start set timer.
        setTic = tic;    

        % Perform analysis
        if strcmpi(correlationType, 'fmc')
            fmiErrorAnalysisMonteCarlo(JOBFILE, SAVEPATH, imageFilePath, imageDigits, parametersPath, SPATIALWINDOW, IMAGESPECTRALFILTER, consts, nProcessors); 
        elseif strcmpi(correlationType, 'rpc')
            rpcErrorAnalysisMonteCarlo(JOBFILE, IMDIR, SAVEPATH, imageFilePath, imageDigits, parametersPath, SPATIALWINDOW, IMAGESPECTRALFILTER, consts, SAVEDATA, nProcessors); 
        elseif strcmpi(correlationType, 'lpc');
            lpcErrorAnalysisMonteCarlo(JOBFILE, SAVEPATH, imageFilePath, imageDigits, parametersPath, SPATIALWINDOW, IMAGESPECTRALFILTER, consts, nProcessors); 
        elseif strcmpi(correlationType, 'strip')
            fmcStripErrorAnalysisMonteCarlo(JOBFILE, SAVEPATH, imageFilePath, parametersPath, SPATIALWINDOW, consts, nProcessors); 
        end

        % Display elapsed time.
        setTime = toc(setTic);
        fprintf(1, '%0.2f sec\n', setTime);

    else
    disp(['Skipping set ' num2str(setList(k))]); 
    end


end
    
end % End if  

end




    
    

