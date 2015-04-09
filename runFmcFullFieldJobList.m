function runFmcFullFieldJobList(JOBLIST)

% Add compiled code path
addpath ba_interpolation;

% Determine the number of jobs in the job list.
nJobs = length(JOBLIST);

% Run each job in the job list
for n = 1 : nJobs
    
    job_tic = tic;
    
    % Extract the jobfile from the joblist.
    JobFile = JOBLIST(n);
    
    % Determine whether to skip image pairs whose results exist.
    skipExisting = JobFile.JobOptions.SkipExisting;
    
    % Number of processors
    dataRepositoryIsAbsolute = JobFile.JobOptions.DataRepositoryIsAbsolute;
    
    % Read the data repository location.
    dataRepository = JobFile.DataRepository;
    
    % Job information
    imageType = JobFile.ImageType; % 'experimental' or 'synthetic';
    setType = JobFile.SetType; % 'vortexring', 'mc', etc.;
    caseName = JobFile.CaseName; % e.g. 'vortexring_2013-11-12_d03_f60_t06'
    
    % Image parameters
    imageBaseName = JobFile.Parameters.Images.BaseName;
    imageExtension = JobFile.Parameters.Images.Extension;
    numberOfDigits = JobFile.Parameters.Images.NumberOfDigits;
    startImage = JobFile.Parameters.Images.Start;
    endImage = JobFile.Parameters.Images.End;
    frameStep = JobFile.Parameters.Images.FrameStep;
    correlationStep = JobFile.Parameters.Images.CorrelationStep;
    particleConcentration = JobFile.Parameters.Images.ParticleConcentration; % Particle concentration (for synthetic images)

    % Grid parameters
    gridSpacingX = JobFile.Parameters.Processing(1).Grid.Spacing.X;
    gridSpacingY = JobFile.Parameters.Processing(1).Grid.Spacing.Y;
   
    % Interrogation Region Parameters
    regionHeight = JobFile.Parameters.Processing(1).InterrogationRegion.Height;
    regionWidth = JobFile.Parameters.Processing(1).InterrogationRegion.Width; 

    % Correlation parameters
    correlationMethod = JobFile.Parameters.Processing(1).Correlation.Method;

    % Image sets
    startSet = JobFile.Parameters.Sets.Start;
    endSet   = JobFile.Parameters.Sets.End;
    setVector = startSet : endSet;
    nSets = length(setVector);
  
    % Determine if the data is synthetic
    isSynthetic = ~isempty(regexpi(imageType, 'syn'));
    
    % Determine the best fft algoritm to use. 
    fftw('planner', 'patient');
    
    % Loop over the sets.
   
    for s = 1 : nSets;
    
    % Determine the current set number.
    currentSet = setVector(s);
    
    % Specify the directories
    if dataRepositoryIsAbsolute
        imageDirectory = dataRepository;
    else
        % Determine the local path to the project repository
        projectRepository = determineLocalRepositoryPath;

        if isSynthetic
            % Synthetic images, where the particle concentration is specified.
            imageDirectory = fullfile(projectRepository, dataRepository, imageType, setType, caseName, ['c_' num2str(particleConcentration, '%0.4f')], [imageBaseName num2str(currentSet, '%05.0f')], 'raw');
        else
            % Experimental images, where the images are just in a folder called "raw"
            imageDirectory = fullfile(projectRepository, dataRepository, imageType, setType, caseName, 'raw');
        end
    end
    
    % Determine the names of output files and directories
  
    % Determine the path to the output directory and create it if it doesn't exist.
    outputDirectory = fullfile(imageDirectory, '..', 'vect', correlationMethod, [num2str(regionHeight) 'x' num2str(regionWidth)], ['cstep_' num2str(correlationStep, '%02.0f')]);
    if ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory);
    end
    
    plotDirectory = fullfile(imageDirectory, '..', 'plot', correlationMethod, [num2str(regionHeight) 'x' num2str(regionWidth)], ['cstep_' num2str(correlationStep, '%02.0f')]);
    if ~exist(plotDirectory, 'dir')
        mkdir(plotDirectory);
    end
    
    % Determine the output name for the vector fields
    outputFileBaseName = [imageBaseName correlationMethod '_grid' num2str(gridSpacingX) 'x' num2str(gridSpacingY) '_region' num2str(regionHeight) 'x' num2str(regionWidth) '_'];
  
    % String specifying the number format for the images.
    numberFormat = ['%0' num2str(numberOfDigits) '.0f'];
    
% Run the job

    % Build a list of image numbers
    firstImageNumbers = startImage : frameStep : endImage;
    secondImageNumbers = firstImageNumbers + correlationStep;
    
    % Determine the number of images
    nPairs = length(firstImageNumbers);
    
    % Build a list of image file paths
    for k = 1 : nPairs
        
        firstImageFilePaths(k, :) = fullfile(imageDirectory, [imageBaseName num2str(firstImageNumbers(k), numberFormat) imageExtension]);
        secondImageFilePaths(k, :) = fullfile(imageDirectory, [imageBaseName num2str(secondImageNumbers(k), numberFormat) imageExtension]);
        outputFilePath(k, :) = fullfile(outputDirectory, [outputFileBaseName num2str(firstImageNumbers(k), numberFormat) '_' num2str(secondImageNumbers(k), numberFormat) '.mat']);
        
        FilePaths(k).FirstImagePath = firstImageFilePaths(k, :);
        FilePaths(k).SecondImagePath = secondImageFilePaths(k, :);
        FilePaths(k).OutputFilePath = outputFilePath(k, :);
        
    end

    % Do the correlations
    for k = 1 : nPairs
        
        if ~(skipExisting && exist(FilePaths(k).OutputFilePath, 'file'))
            disp(['Correlating pair ' num2str(k) ' of ' num2str(nPairs)]);
            pair_tic = tic;
            fmcFullField(FilePaths(k), JobFile);
            disp(['Saved vector field to ' outputFilePath(k, :)]);
            disp(['Image Pair Time: ' num2str(toc(pair_tic)) ' sec'])
        end
    end
   
    
    end % End looping over sets   

    job_toc = toc(job_tic);
    fprintf('Total job time: %d seconds\n', job_toc);
    
end


end


