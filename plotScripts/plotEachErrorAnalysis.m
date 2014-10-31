% function runFmiJobFileMc(SETTYPE, CASENAME, DIMENSIONS, STARTSET, ENDSET, NPROCESSORS)
function plotEachErrorAnalysis(JOBLIST)

nJobs = length(JOBLIST);

for n = 1 : nJobs
    
JOBFILE = JOBLIST(n);

SAVEDATA = 1; % Force data save (hard code so that I don't spend a week processing data and forget to save the results)

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

% Main code repository
Repository = determineLocalRepositoryPath;

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

% Close any open files
fclose all;

% Create error log file
fid = fopen(fullfile( '~/Desktop', 'FMIerrorLog.txt' ), 'w');

% Print header to error log file
fprintf( fid, [ 'FMI-RPC Processing Error Log \nCase: ' setBase '\n\n' ] ); 
fclose(fid);

    for k = startSet :  endSet
        
    fprintf(1, ['Analyzing set ' ...
             setBase num2str(k, setFormat) ' (' num2str(k - startSet + 1)...
             ' of ' num2str(endSet - startSet + 1) ')\n' ], k); % Display status
         
         
         SAVEPATH = fullfile( writeDir, [ saveBase num2str( k, setFormat ) '.mat' ] ); % Save path
         
        load(SAVEPATH);
        
        rotDeg = rad2deg(trueRotation);
        
        translationErrorMagnitude = sqrt(translationAbsErrorX.^2 + translationAbsErrorY.^2);
        
        plot(rotDeg, translationErrorMagnitude, 'ok'); 
        xlim([0 360]);
        ylim([0 20]);
        drawnow;
        if max(translationErrorMagnitude(:)) > 10
            keyboard
        end
        hold on;
        
       
    end
    
end % End if  

end




    
    

