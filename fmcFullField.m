
function fmcFullField(FILEPATHS, JOBFILE)

% Check for whether to use compiled codes
% Only use compiled codes on Linux.
% This can be changed if you want to re-compile the codes
% on your own platform, i.e. if you figure out
% how to compile mex files on Mac, etc.
COMPILED = JOBFILE.JobOptions.RunCompiled;

% Number of PIV passes
% Take the minimum of the requested number of passes and the number of pass
% parameter structures specified.
numberOfPasses = min(JOBFILE.JobOptions.NumberOfPasses, length(JOBFILE.Parameters.Processing));

% Chop the jobfile processing parameters down to the first n passes
JOBFILE.Parameters.Processing = JOBFILE.Parameters.Processing(1 : numberOfPasses);

% Read in the clean images.
image1_raw = double(imread(FILEPATHS.FirstImagePath));
image2_raw = double(imread(FILEPATHS.SecondImagePath));

% Image dimensions
[imageHeight, imageWidth] = size(image1_raw);

% Check whether to re-start from a previously existing velocity field.
startFromExistingField = JOBFILE.JobOptions.StartFromExistingField;

% If starting from a previously existing velocity field,
% load that velocity field.
if startFromExistingField && exist(FILEPATHS.OutputFilePath, 'file') 
    
    % Load the existing velocity field if it exists.
    load(FILEPATHS.OutputFilePath);
        
    % Set the current pass number to the first specified pass
    % If "start from existing" path is specified and 
    % a starting-pass number is specified less than 2,
    % then default to starting at pass 2. This 
    % deals with inputs of 0, negative numbers, 
    % or confusion about what "starting pass" means.
    % Presumably if the user specifies starting from an existing
    % pass, they want to pick things up starting with pass 2
    % (starting on pass 1 would just repeat pass 1, which is the same
    % thing as startFromExistingField = 0)
    if JOBFILE.JobOptions.StartPass < 2
        thisPass = 2;
    else
        thisPass = JOBFILE.JobOptions.StartPass;
    end
    
    % This sets a loop counter for the PIV passes.
    p = thisPass - 1;
    
    % Create the function variables from the saved data
    % This just entails flipping all the coordinate and 
    % displacement data
    for pass = 1 : length(X)
        gx{pass} = flipud(X{pass});
        gy{pass} = flipud(Y{pass});
        TRANSLATIONX{pass} = flipud(U{pass});
        TRANSLATIONY{pass} = flipud(V{pass});
        ROTATION{pass} = flipud(R{pass});
        SCALING{pass} = flipud(S{pass});
        uVal{pass} = flipud(UVAL{pass});
        vVal{pass} = flipud(VVAL{pass});
    end
    
else
    % Initialize the counter for the number of user-specified passes
    thisPass = 1;
    p = 0;
end

% Overwrite the previous pass jobfile and filename info.
FilePaths = FILEPATHS;
JobFile = JOBFILE;

% Initialize the DWO iteration counter (this is for DWO convergence)
nDwoIterations = 0;

% Loop over the passes.
% for p = 1 : numberOfPasses;
while thisPass <= numberOfPasses;
    
    % Increment the pass counter 
    p = p + 1;
    
    % Write the number of the specified pass that's currently executing
    PASSNUMBER(p) = thisPass;
    
    % Read deformation flag.
    doImageDeformation = JobFile.Parameters.Processing(p).DoImageDeformation;
    
    % Read smoothing parameters
    doSmoothing = JobFile.Parameters.Processing(p).Smoothing.DoSmoothing;
    
    % Smoothing kernel diameter
    smoothingKernelDiameter = JobFile.Parameters.Processing(p).Smoothing.KernelDiameter;
   
    % Smoothing kernel gaussian standard deviation  
    smoothingGaussianStdDev = JobFile.Parameters.Processing(p).Smoothing.KernelGaussianStdDev;
    
    % Flag to calculate image disparity
    calculateImageDisparity = JobFile.Parameters.Processing(p).DoImageDisparity;
    
    % Check deformation flag
    if p > 1 && doImageDeformation
        
        % Smooth field if specified
        if doSmoothing
            % Smooth the velocity field.
            source_field_u{p-1} = smoothField(uVal{p-1}, smoothingKernelDiameter, smoothingGaussianStdDev);
            source_field_v{p-1} = smoothField(vVal{p-1}, smoothingKernelDiameter, smoothingGaussianStdDev);
       
        else
            source_field_u{p-1} = uVal{p-1};
            source_field_v{p-1} = vVal{p-1};
        end
        
        % Create the pixel coordinates.
        [xi_integer, yi_integer] = meshgrid(1:imageWidth, 1:imageHeight);
        
        % Shift the pixel coordinates by 0.5 pixels
        XI = xi_integer - 0.5;
        YI = yi_integer - 0.5;
        
        % Create interpolation structures for the velocity field.
        % Temporary: change from spline to cubic interpolation, and change
        % from lienar to nearest neighbor extrapolation. This is for
        % comparison with prana.
        interpolant_tx = griddedInterpolant(gy{p-1}, gx{p-1}, source_field_u{p-1}, 'cubic', 'nearest');
        interpolant_ty = griddedInterpolant(gy{p-1}, gx{p-1}, source_field_v{p-1}, 'cubic', 'nearest');
        
        % This is the velocity field upsampled to every pixel.
        UI = interpolant_tx(YI, XI);
        VI = interpolant_ty(YI, XI);
        
        % These are the coordinates at which to resample image 1.
        XD1 = XI - UI/2;
        YD1 = YI - VI/2;
        
        % These are the coordinates at which to resample image 2.
        XD2 = XI + UI/2;
        YD2 = YI + VI/2;

        % Resample the images
        image1 = sincBlackmanInterp2(image1_raw, XD1 + 0.5, YD1 + 0.5, 8, 'blackman');
        image2 = sincBlackmanInterp2(image2_raw, XD2 + 0.5, YD2 + 0.5, 8, 'blackman');
        
    else
        
        % If not deform or if we're on the first pass, use the raw images.
        image1 = image1_raw;
        image2 = image2_raw;
    end
    
    % Interrogation region dimensions
    regionHeight = JobFile.Parameters.Processing(p).InterrogationRegion.Height;
    regionWidth = JobFile.Parameters.Processing(p).InterrogationRegion.Width;    
    
    % Fmc difference method
    fmcDifferenceMethod_string = JobFile.Parameters.Processing(p).FmcDifferenceMethod; 
    
    % Determine what FMC difference method to use. 
    % fmcDifferenceMethod = 1 is central difference.
    % fmcDifferenceMethod = 2 is forward difference.
    % fmcDifferenceMethod = 3 is backward difference.
    if ~isempty(regexpi(fmcDifferenceMethod_string, 'bac'))
        fmcDifferenceMethod = 3;
    elseif ~isempty(regexpi(fmcDifferenceMethod_string, 'for'))
        fmcDifferenceMethod = 2;
    else
        % Default to central difference
        fmcDifferenceMethod = 1;
    end
    
    % Apodization window parameters.
    spatialWindowFraction = JobFile.Parameters.Processing(p).InterrogationRegion.SpatialWindowFraction;
    fmiWindowSize = JobFile.Parameters.Processing(p).InterrogationRegion.FMIWindowSize;
    fmiWindowType = JobFile.Parameters.Processing(p).InterrogationRegion.FMIWindowType;

    % FFT Parameters
    fftSize = JobFile.Parameters.Processing(p).FFTSize;
    spectrum_height = fftSize(1);
    spectrum_width  = fftSize(2);

    % Image resampling parameters
    numberOfRings = JobFile.Parameters.Processing(p).Resampling.NumberOfRings;
    numberOfWedges = JobFile.Parameters.Processing(p).Resampling.NumberOfWedges;
    rMin = JobFile.Parameters.Processing(p).Resampling.MinimumRadius;
    rMax = min(spectrum_height, spectrum_width) / 2 - 1;

    % Correlation parameters
    correlationMethod = JobFile.Parameters.Processing(p).Correlation.Method;
    
    % Determine which correlation type to use. Fmc or RPC
    isFmc = ~isempty(regexpi(correlationMethod, 'fmc'));
    isRpc = ~isempty(regexpi(correlationMethod, 'rpc'));
    isScc = ~isempty(regexpi(correlationMethod, 'scc'));
    
    % RPC diameters
    spatialRPCDiameter = JobFile.Parameters.Processing(p).Correlation.SpatialRPCDiameter;
    fmiRpcDiameter = JobFile.Parameters.Processing(p).Correlation.FMCDiameter; 

    % Create the gaussian intensity window to be applied to the the raw image interrogation regions
    spatialWindow = gaussianWindowFilter_prana([regionHeight, regionWidth], spatialWindowFraction .* [regionHeight, regionWidth]);
    
    % Determine the FMI window type.
    isHann1 = ~isempty(regexpi(fmiWindowType, 'hann1'));
    isHann2 = ~isempty(regexpi(fmiWindowType, 'hann2'));
    isGaussianSkew = ~isempty(regexpi(fmiWindowType, 'gauss_skew'));

    % Create the FMI Window
    if isHann1
        fmiWindow1D = hann1(numberOfRings, [fmiWindowSize(1) fmiWindowSize(2)], fmiWindowSize(3));
        fmiWindow = repmat(fmiWindow1D, numberOfWedges, 1);
    elseif isHann2
        fmiWindow = hann2([numberOfWedges, numberOfRings], fmiWindowSize(1));
    elseif isGaussianSkew
        fmiWindow1D = gaussianWindowFilter_asymmetric(numberOfRings, fmiWindowFraction);
        fmiWindow = repmat(fmiWindow1D, numberOfWedges, 1);
    else
        fmiWindow = gaussianWindowFilter([numberOfWedges, numberOfRings], fmiWindowSize, 'fraction');
    end

    % Create the gaussian spectral energy filter be applied to the raw image correlation
    imageSpectralFilter = spectralEnergyFilter(regionHeight, regionWidth, spatialRPCDiameter); 

    % Create the FMI spectral filter (i.e. the FMI RPC filter).
    fmiSpectralFilter = spectralEnergyFilter(numberOfWedges, numberOfRings, fmiRpcDiameter); 

    % Make a matrix of the subregion coordinates.
    % Do this only once to increase speed (meshgrid is slow).
    [xImage, yImage] = meshgrid(1 : regionWidth, 1 : regionHeight);

    % FFT Spectrum coordinates
    [xSpectrum, ySpectrum] = meshgrid(1:spectrum_width, 1:spectrum_height);

    % Figure out the log polar resampling coordinates.
    [xLP, yLP] = LogPolarCoordinates([spectrum_height, spectrum_width], numberOfWedges, numberOfRings, rMin, rMax, 2 * pi);

    % Universal Outlier Detection Parameters
    uodStencilRadius = JobFile.Parameters.Processing(p).Validation.UodStencilRadius;
    uodThreshold = JobFile.Parameters.Processing(p).Validation.UodThreshold;
    uodExpectedDifference = JobFile.Parameters.Processing(p).Validation.UodExpectedDifference;
    
    % Save the image size to the processing field for easy passing around.
    JobFile.Parameters.Processing(p).Images.Height = imageHeight;
    JobFile.Parameters.Processing(p).Images.Width  = imageWidth;

    % Flag specifying whether or not to do DWO
    doDiscreteWindowOffset = JobFile.Parameters.Processing(p).DoDiscreteWindowOffset;
    
    % Discrete window offset differencing method
    dwoDifferenceMethod = JobFile.Parameters.Processing(p).DwoDifferenceMethod;
    
    % Grid parameters
    gridSpacingX = JobFile.Parameters.Processing(p).Grid.Spacing.X;
    gridSpacingY = JobFile.Parameters.Processing(p).Grid.Spacing.Y;

    % Make sure the grid buffer is at least half the size of the interrogation region
    gridBufferY = max(JobFile.Parameters.Processing(p).Grid.Buffer.Y, ceil(regionHeight / 2));
    gridBufferX = max(JobFile.Parameters.Processing(p).Grid.Buffer.X, ceil(regionWidth / 2));

    % If the pass number is greater than one, i.e., if at least one pass has
    % finished, and also if discrete window offset is enabled,
    % then interpolate the velocity field from the previous pass
    % onto the grid for the current pass.
    % Round the grid shift values so that grid points are shifted
    % from integer coordinates to integercoordinates
    if p > 1 && doDiscreteWindowOffset
        [gx{p}, gy{p}, gx_01, gy_01, gx_02, gy_02] = discreteWindowOffset(gx{p-1}, gy{p-1}, uVal{p-1}, vVal{p-1}, JobFile.Parameters.Processing(p));
        
    else
         % Generate the list of coordinates that specifies the (X, Y) centers of all of the interrogation regions 
        [ gx{p}, gy{p} ] = gridImage([imageHeight, imageWidth], [gridSpacingY gridSpacingX], gridBufferY, gridBufferX);
   
        % If this is the first pass or if DWO is not specified, then keep
        % the original grid points.
        gx_01 = gx{p};
        gy_01 = gy{p};
        gx_02 = gx{p};
        gy_02 = gy{p}; 
    end
    
    % Calculate the grid shifts for both images. 
    % These should all be integers!!
    gridShiftX_01 = gx_01 - gx{p};
    gridShiftY_01 = gy_01 - gy{p};
    gridShiftX_02 = gx_02 - gx{p};
    gridShiftY_02 = gy_02 - gy{p};

    % Determine the size of the grid (number of grid points)..
    [numRows, numColumns] = size(gx_01);

    % Determine the number of interrogation regions to be correlated
    nRegions = numRows * numColumns;
    
    %%%% TEMP
    grid_x_temp_01{p} = gx_01(:);
    grid_y_temp_01{p} = gy_01(:);

    grid_x_temp_02{p} = gx_02(:);
    grid_y_temp_02{p} = gy_02(:);

    % Extract the subregions from each image.
    regionMatrix1 = extractSubRegions(image1, [regionHeight, regionWidth], gx_01(:), gy_01(:));
    regionMatrix2 = extractSubRegions(image2, [regionHeight, regionWidth], gx_02(:), gy_02(:));
    
    % Preallocate memory for the vectors to hold the estimates of translation, rotation, and scaling.
    estimatedTranslationY = zeros(nRegions, 1);
    estimatedTranslationX = zeros(nRegions, 1);
    estimatedRotation = zeros(nRegions, 1); 
    estimatedScaling = ones(nRegions, 1);
    fmiTranslationY = zeros(nRegions, 1);
    fmiTranslationX = zeros(nRegions, 1);
    
    % Preallocate memory for disparity
    disparity_x_vector = zeros(nRegions, 1);
    disparity_y_vector = zeros(nRegions, 1);
    
    % Number of particles identified in each region
    % this is for image disparity calculation
    nParticles = zeros(nRegions, 1);

    % Initialize FMC peak height ratio vector.
    fmcPeakRatio = zeros(nRegions, 1);

    % Initialize RPC peak height ratio vector.
    spatialPeakRatio = zeros(nRegions, 1);
    
    % height of the primary spatial peaks
    spatialPeakHeight = zeros(nRegions, 1);
    
    % Diameter of primary spatial peaks
    spatialPeakDiameter = zeros(nRegions, 1);
    
    % Start a timer
    t = tic;
    
    % Do all the correlations for the image.
    parfor k = 1 : nRegions
        
        % Extract the subregions from the subregion stacks.
        subRegion1 = regionMatrix1(:, :, k);
        subRegion2 = regionMatrix2(:, :, k);
                
        % Perform FMC processing. 
        if isFmc
            % Perform the FMC correlation.
            [estimatedTranslationY(k), estimatedTranslationX(k),...
            estimatedRotation(k), estimatedScaling(k), ...
            fmcPeakRatio(k), spatialPeakRatio(k), spatialPeakHeight(k), spatialPeakDiameter(k)] = ...
            ...
            FMC(subRegion1, subRegion2, spatialWindow, imageSpectralFilter,...
            fmiWindow, fmiSpectralFilter, ...
            xImage, yImage, ...
            spectrum_height, spectrum_width,...
            xLP, yLP, rMin, rMax, fmcDifferenceMethod, COMPILED);
        
        % Perform RPC analysis 
        % The zero in this input means "Do not search multiple peaks,"
        % i.e., use only the primary peak.
        elseif isRpc
            [estimatedTranslationY(k), estimatedTranslationX(k), rpcPlane, spatialPeakHeight(k), spatialPeakDiameter(k)]...
                = RPC(spatialWindow .* subRegion1, spatialWindow .* subRegion2,...
                imageSpectralFilter, COMPILED); 

            % Measure the peak height ratio
            if COMPILED
                spatialPeakRatio(k) = measurePeakHeightRatio(rpcPlane, COMPILED);
            else
                spatialPeakRatio(k) = measurePeakHeightRatio(rpcPlane, COMPILED);
            end

        % Perform SCC analysis.
        elseif isScc
            [estimatedTranslationY(k), estimatedTranslationX(k), spatialPeakRatio(k)]...
                = SCC(spatialWindow .* subRegion1, spatialWindow .* subRegion2);
        end
        
        % These lines will calculate the disparity between regions.
        if calculateImageDisparity
            [DISPARITY_X, DISPARITY_Y] = calculateTransformedRegionDisparity(subRegion1, subRegion2,...
                estimatedTranslationY(k), estimatedTranslationX(k),...
                estimatedRotation(k), estimatedScaling(k), xImage, yImage, 0.95, 2, COMPILED);
        end
        
    end % end for k = 1 : nRegions

    % Inform the user
    disp(['Correlation times (pass ' num2str(p) '): ' num2str(toc(t)) ' sec' ]);
    disp('');
    
    % Reshape the raw measured displacements into matrices.
    tx_raw{p} = reshape(estimatedTranslationX, numRows, numColumns);
    ty_raw{p} = reshape(estimatedTranslationY, numRows, numColumns);
    
    % Shift the measured velocities by the deform or DWO values. 
    if doImageDeformation
        
        % Only shift for the second and greater iterations
        if p > 1
            
            % Temporary change to cubic and nearest
            % Create interpolant structures
              interpolant_tx = griddedInterpolant(YI, XI, UI, 'cubic', 'nearest');
              interpolant_ty = griddedInterpolant(YI, XI, VI, 'cubic', 'nearest');

            % Evaluate the interpolant structures
              tx_shift = interpolant_tx(gy{p}, gx{p});
              ty_shift = interpolant_ty(gy{p}, gx{p});

        else
            % For the first pass, set the shift values to zero.
            tx_shift = zeros(size(gx{p}));
            ty_shift = zeros(size(gx{p}));
        end
        
        % Add the shift values to the measured displacements for deform.
        TRANSLATIONX{p} = tx_raw{p} + tx_shift;
        TRANSLATIONY{p} = ty_raw{p} + ty_shift;
        
    else     
        % If not deform then just shift the velocities by the DWO grid
        % shift values, which are already zeros if DWO wasn't used.
        TRANSLATIONX{p} = tx_raw{p} + gridShiftX_02 - gridShiftX_01;
        TRANSLATIONY{p} = ty_raw{p} + gridShiftY_02 - gridShiftY_01;
    end
    
    % Reshape the rotation and scaling measurements into matrices.
    ROTATION{p} = reshape(estimatedRotation, numRows, numColumns);
    SCALING{p} =  reshape(estimatedScaling, numRows, numColumns);

    % Reshape the peak ratio measurements into matrices.
    SPATIAL_PEAK_RATIO{p} = flipud(reshape(spatialPeakRatio, numRows, numColumns));
    FMC_PEAK_RATIO{p} = flipud(reshape(fmcPeakRatio, numRows, numColumns));
    
    % Reshape the disparity vectors into matrices.
    DISPARITY_X{p} = flipud(reshape(disparity_x_vector, numRows, numColumns));
    DISPARITY_Y{p} = flipud(reshape(disparity_y_vector, numRows, numColumns));
    N_PARTICLES{p} = flipud(reshape(nParticles, numRows, numColumns));

    % Run Prana's validation code. Note that right now the rotation
    % The previous codes to do this were "universalOutlierDetection.m"
    % and "universalOutlierReplacement.m"
    % estimate isn't validated.
    [uVal{p}, vVal{p}, isOutlier{p}] = validateField_prana(gx{p}, gy{p}, TRANSLATIONX{p}, TRANSLATIONY{p}, uodExpectedDifference);
        
    % Check for convergence if it's requested. 
    % If the velocity estimate has converged, go on
    % to the next user-specified pass. Otherwise, 
    % repeat the previous pass. 
    
    % This determines whether DWO convergence iterations were specified
    doDwoConvergence = JobFile.Parameters.Processing(p).DwoConverge;
    
    % If DWO convergence was specified (and at least one pass has
    % completed), then check the other parameters regarding DWO convergence
    if doDwoConvergence
        % Inform the user
        disp('Convergence requested. Checking convergence.')
        
        % Determine the maximum number of iterations specified
        maxDwoConvergenceIterations = JobFile.Parameters.Processing(p).DwoMaxConvergenceIterations;
        
        % Determine the DWO convergence criteria
        dwoConvergenceCriteria = JobFile.Parameters.Processing(p).DwoConvergenceCriteria;
        
        % If at least one DWO iteration has been completed...
        if nDwoIterations > 0
            % Determine the 2-norm of the velocity field components compared to
            % the previous pass. This is the metric against which the
            % convergence criteria is compared.
            uNorm(p) = mean(abs(TRANSLATIONX{p}(:) - TRANSLATIONX{p-1}(:)));
            vNorm(p) = mean(abs(TRANSLATIONY{p}(:) - TRANSLATIONY{p-1}(:)));
            rNorm(p) = mean(abs(ROTATION{p}(:) - ROTATION{p-1}(:)));
            
            % Inform the user
            disp(['U norm: ' num2str(uNorm(p), '%10.3e') '    V norm: ' num2str(vNorm(p), '%10.3e') '    Criteria: ' num2str(dwoConvergenceCriteria, '%10.3e')]);

            % Check if the convergence criteria have been reached
            hasConverged(p) = min(uNorm(p) <= dwoConvergenceCriteria , vNorm(p) <= dwoConvergenceCriteria);
        else
            % Convergence is never reached before the first DWO iteration.
            hasConverged(p) = 0;
        end
        
        % If the velocity field has converged or the max number of
        % iterations has been reached
        if hasConverged(p)
           
            % Inform the user
            disp(['Pass ' num2str(thisPass) ' converged after ' num2str(nDwoIterations) ' DWO iterations. Incrementing pass.']);
            
            % Save the number of iterations
            dwoIterations(thisPass) = nDwoIterations;
     
            % Reset the DWO Iteration counter
            nDwoIterations = 0;
            
            % Increment the counter for the user-specified passes
            thisPass = thisPass + 1;
           
        % Max iterations reached.
        elseif nDwoIterations > maxDwoConvergenceIterations
            % Inform the user
            disp(['Max number of iterations reached for ' num2str(thisPass) '.']);
            
            % Save the number of iterations
            dwoIterations(thisPass) = nDwoIterations;
     
            % Reset the DWO Iteration counter
            nDwoIterations = 0;
            
            % Increment the counter for the user-specified passes
            thisPass = thisPass + 1;
        
        % Neither convergence nor max iterations have been reached.    
        else 
            
            % Inform the user
            disp(['Pass ' num2str(thisPass) ' has not converged after ' num2str(nDwoIterations)  ' DWO iterations. Iterating DWO.']);
            
            % Increment the DWO iteration counter
            nDwoIterations = nDwoIterations + 1; 

            % Update the jobfile with a new pass
            % Then exit the convergence-checking loop with the updated jobfile
            JobFile.Parameters.Processing(p + 1 : end + 1) = JobFile.Parameters.Processing(p : end);
            
        end
        
    else
        % Increment the counter for the user-specified passes
        thisPass = thisPass + 1;
        hasConverged(p) = 0;
    end
        
end % End for p = 1 : numberOfPasses  

% Number of passes that ended up getting run
finalNumberOfPasses = length(JobFile.Parameters.Processing);

% This flips everything over the Y axis and saves the variables to the
% output structures.
for p = 1 : finalNumberOfPasses
    X{p} = flipud(gx{p});
    Y{p} = flipud(gy{p});
    U{p} = flipud(TRANSLATIONX{p});
    V{p} = flipud(TRANSLATIONY{p});
    R{p} = flipud(ROTATION{p});
    S{p} = flipud(SCALING{p});
    UVAL{p} = flipud(uVal{p});
    VVAL{p} = flipud(vVal{p});
    IS_OUTLIER{p} = flipud(isOutlier{p});
    
    % TEMPORARY: Don't validate rotation estimate.
    RVAL{p} = zeros(size(R{p}));
    
end

% Rename variable
CONVERGED = hasConverged;

% Set source field variables to zeros if only one pass was specifed.
if numberOfPasses < 2
    source_field_u{1} = zeros(size(X{1}));
    source_field_v{1} = zeros(size(X{1}));
end

% Save the results
save(FilePaths.OutputFilePath, ...
    'X', 'Y', 'U', 'V', 'R', 'S', 'IS_OUTLIER',...
    'UVAL', 'VVAL', 'RVAL', 'tx_raw', 'ty_raw', 'DISPARITY_X', 'DISPARITY_Y', 'N_PARTICLES',...
    'FMC_PEAK_RATIO', 'SPATIAL_PEAK_RATIO', 'PASSNUMBER', 'CONVERGED', 'source_field_u', 'source_field_v', ...
    'FilePaths', 'JobFile');

end






















