function JobInfo = determineJobFileNames(JOBLIST)
% This function returns all the file names / directories 
% for a job list, including raw images, vector files, and error files.


%% Determine the local path to the project repository
projectRepository = determineLocalRepositoryPath;

% Determine the number of jobs in the job list.
nJobs = length(JOBLIST);

% Initialize all the job info structures
JobInfo(nJobs).Sets = [];

% for n = 1 : nJobs
%     startSet = JOBLIST(n).Parameters.Sets.Start;
%     endSet   = JOBLIST(n).Parameters.Sets.End;
%     skipSet  = JOBLIST(n).Parameters.Sets.Skip;
%     
%     % Number of sets in the job
%     nSets = length(startSet : skipSet : endSet);
%     JobInfo(n).Sets(nSets) =
%     
%     for s = 1 : nSets
%         
%     end
% end

disp('Determining file names...');

%% Run each job in the job list
for n = 1 : nJobs
    
    %% Parse the jobfile.
    
    % Extract the jobfile from the joblist.
    JobFile = JOBLIST(n);
   
    % Save job information to output structure
    JobInfo(n).JobFile = JobFile;
   
    % Determine if data path is absolute.
    dataRepositoryIsAbsolute = JobFile.JobOptions.DataRepositoryIsAbsolute;
    
    % Job information
    imageType = JobFile.ImageType; % 'experimental' or 'synthetic';
    setType = JobFile.SetType; % 'vortexring', 'mc', etc.;
    caseName = JobFile.CaseName; % e.g. 'vortexring_2013-11-12_d03_f60_t06'
    dataRepository = JobFile.DataRepository;
    
    % Set information
    startSet = JobFile.Parameters.Sets.Start;
    endSet  = JobFile.Parameters.Sets.End;
    skipSet = JobFile.Parameters.Sets.Skip;
    
    % Image parameters
    imageBaseName = JobFile.Parameters.Images.BaseName;
    imageExtension = JobFile.Parameters.Images.Extension;
    numberOfDigits = JobFile.Parameters.Images.NumberOfDigits;
    startImage = JobFile.Parameters.Images.Start;
    endImage = JobFile.Parameters.Images.End;
    frameStep = JobFile.Parameters.Images.FrameStep;
    correlationStep = JobFile.Parameters.Images.CorrelationStep;

    % Grid parameters
    gridSpacingX = JobFile.Parameters.Processing(1).Grid.Spacing.X;
    gridSpacingY = JobFile.Parameters.Processing(1).Grid.Spacing.Y;
%     gridBufferX = JobFile.Parameters.Processing.Grid.Buffer.X;
%     gridBufferY = JobFile.Parameters.Processing.Grid.Buffer.Y;
   
    % Interrogation Region Parameters
    regionHeight = JobFile.Parameters.Processing(1).InterrogationRegion.Height;
    regionWidth = JobFile.Parameters.Processing(1).InterrogationRegion.Width; 
%     spatialWindowSize = JobFile.Parameters.Processing(1).InterrogationRegion.SpatialWindowDimensions;
%     spatialWindowType = JobFile.Parameters.Processing(1).InterrogationRegion.SpatialWindowType;
%     fmiWindowSize = JobFile.Parameters.Processing.InterrogationRegion.FMIWindowDimensions;
%     fmiWindowType = JobFile.Parameters.Processing.InterrogationRegion.FMIWindowType;
%     
    % List of sets
    setList = startSet : skipSet : endSet;
    
    % Number of sets
    nSets = length(setList);
    
    % Particle concentration (for synthetic images)
    particleConcentration = JobFile.Parameters.Images.ParticleConcentration;
    
    % Correlation parameters
    correlationMethod = JobFile.Parameters.Processing(1).Correlation.Method;

    % String specifying the number format for the images.
    numberFormat = ['%0' num2str(numberOfDigits) '.0f'];
    
    % List of the image numbers used in each correlation pair.
    firstImageNumbers = startImage : frameStep : endImage;
    secondImageNumbers = firstImageNumbers + correlationStep;
    
    % Number of image pairs
    nPairs = length(firstImageNumbers);
   
    % Determine the output name for the vector fields
    JobInfo(n).VectorFileBaseName = [imageBaseName correlationMethod '_grid' num2str(gridSpacingX) 'x' num2str(gridSpacingY) '_region' num2str(regionHeight) 'x' num2str(regionWidth) '_'];
  
    % Determine the base name for the error plots
    JobInfo(n).ErrorFileBaseName = [imageBaseName correlationMethod 'error_grid' num2str(gridSpacingX) 'x' num2str(gridSpacingY) '_region' num2str(regionHeight) 'x' num2str(regionWidth) '_'];
    
    %% Determine file paths

    % Flag the images as 'synthetic' if they are specified as such.
    isSynthetic = ~isempty(regexpi(imageType, 'syn'));
    
    % Determine directory paths for synthetic images.
    if dataRepositoryIsAbsolute
        caseDirectory = dataRepository;
    else
        if isSynthetic % Synthetic images, relative repository path
            caseDirectory = fullfile(projectRepository, dataRepository, imageType, setType, caseName, ['c_' num2str(particleConcentration, '%0.4f')]);
        else % Experimental images, relative repository path
            caseDirectory = fullfile(projectRepository, dataRepository, imageType, setType, caseName);
        end
    end
    
    % Statistics directory
    JobInfo(n).StatsDir = fullfile(caseDirectory, 'stat');
    
    % Loop over all the image sets
    for s = 1 : nSets
        
        % Determine the number of the current set
        currentSet = setList(s);
        
        JobInfo(n).Sets(s).SetNumber = currentSet;
        
        % Path to the set directory   
        setDirectory = fullfile(caseDirectory, [imageBaseName num2str(currentSet, '%05.0f')]);
        
        % Imate directory
        imageDirectory = fullfile(setDirectory, 'raw');

        % Path to the directory containing the simulation parameters
        parametersDirectory = fullfile(setDirectory, 'parameters');

        % Path to the file containing the parameters
        parametersPath = fullfile(parametersDirectory, 'parameters.mat');

        % Path to the directory containing the PIV output vectors.
        vectorDirectory = fullfile(setDirectory, 'vect', correlationMethod, [num2str(regionHeight) 'x' num2str(regionWidth)], ['cstep_' num2str(correlationStep, '%02.0f')]);

        % Path to the directory to which to write the error analysis data.
        % Make this if it doesn't already exist.
        errorDirectory = fullfile(setDirectory, 'err', correlationMethod, [num2str(regionHeight) 'x' num2str(regionWidth)], ['cstep_' num2str(correlationStep, '%02.0f')]);

        % Save names to the output structure.
        JobInfo(n).Sets(s).SetDirectory = setDirectory;
        JobInfo(n).Sets(s).Images.Directory = imageDirectory;
        JobInfo(n).Sets(s).Parameters.Directory = parametersDirectory;
        JobInfo(n).Sets(s).Parameters.Path = parametersPath;
        JobInfo(n).Sets(s).Error.Directory = errorDirectory;
        JobInfo(n).Sets(s).Vector.Directory = vectorDirectory;
    
    end % End for s = 1 : nSets
    

end

    


end


