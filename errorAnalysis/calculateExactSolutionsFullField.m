function calculateExactSolutionsFullField(JOBLIST)

%% Parse the input variables and set relevant variables.
% isDirect = ~isempty(regexpi(VORTMETHOD, 'dir'));
% isGradient = ~isempty(regexpi(VORTMETHOD, 'grad'));
% noGradient = ~(isDirect || isGradient);
% 
% if isDirect
%     gradientMethod = 'Direct';
% elseif isGradient
%     gradientMethod = 'Gradient';
% else
%     gradientMethod = 'Raw';
% end

    
% Vorticity threshold for defining the core-region of the vortex. 
% Regions of the flow where the vorticity magnitude is greater than or equal to 
% this fraction of the maximum vorticity magnitude are defined as being "in the core."
vorticityThreshold = 0.3;

%% Determine the local path to the project repository
projectRepository = determineLocalRepositoryPath;

% Determine the number of jobs in the job list.
nJobs = length(JOBLIST);

% Determine the names for all the files
JobInfo = determineJobFileNames(JOBLIST);



%% Run each job in the job list
for n = 1 : nJobs
    
    %% Parse the jobfile.
    
    % Extract the jobfile from the joblist.
    JobFile = JOBLIST(n);
    
    % Determine the number of passes
    nPasses = length(JobFile.Parameters.Processing);
    
    % Option specifying whether to compare the estimated velocity field
    % to the Eulerian or Lagrangian true velocity field.
    comparisonType = JobFile.JobOptions.ComparisonType;
    isLagrangian = ~isempty(regexpi(comparisonType, 'lag'));
    isEulerian = ~isempty(regexpi(comparisonType, 'eul'));
    
    % Number of sets in this Job
    nSets = length(JobInfo(n).Sets);
    
    % Stats Directory
    statsDir = JobInfo(n).StatsDir;
    if ~exist(statsDir, 'dir');
        mkdir(statsDir)
    end
    
    % Image parameters
    numberOfDigits = JobFile.Parameters.Images.NumberOfDigits;
    startImage = JobFile.Parameters.Images.Start;
    endImage = JobFile.Parameters.Images.End;
    frameStep = JobFile.Parameters.Images.FrameStep;
    correlationStep = JobFile.Parameters.Images.CorrelationStep;
    
    % Determine if the rotation angle was calculated in a left handed
    % coordinate system (i.e. if positive rotation is clockwise).
    isLeftHanded = JobFile.JobOptions.LeftHanded;
    
    % Correlation parameters
    correlationMethod = JobFile.Parameters.Processing(1).Correlation.Method;
    
    % Determine if this is a deform pass
    isDeform = ~isempty(regexpi(correlationMethod, 'def'));
   
    % List of the image numbers used in each correlation pair.
    firstImageNumbers = startImage : frameStep : endImage;
    secondImageNumbers = firstImageNumbers + correlationStep;
   
    % Number of pairs per set.
    nPairs = length(firstImageNumbers);
    
    % Number format
    numberFormat = ['%0' num2str(numberOfDigits) '.0f'];
    
    % Vector file base name
    vectorFileBaseName = JobInfo(n).VectorFileBaseName;
    
    % Error file base name
    errorFileBaseName = JobInfo(n).ErrorFileBaseName;
    
    % Vector directory
    firstVectorDir = JobInfo(n).Sets(1).Vector.Directory;
        
    % Path to the first vector file in the set
    firstVectorName = [vectorFileBaseName num2str(firstImageNumbers(1), numberFormat) '_' num2str(secondImageNumbers(1), numberFormat) '.mat'];
    firstVectorPath = fullfile(firstVectorDir, firstVectorName);
    
    % Load the first vector path.
    firstVectorFile = load(firstVectorPath);
    
    % Average errors inside core
    meanErrInside = zeros(nPairs * nSets, nPasses);
    
    % Calculate the true solutions at the grid points for each pass.  
    for PASS = 1 : nPasses

        % Read the grid points at which vectors were calculated
        X{PASS} = firstVectorFile.X{PASS};
        Y{PASS} = firstVectorFile.Y{PASS};
        
        % Size of the grid
        [gridHeight, gridWidth] = size(X{PASS});  

        % Initialize cumultive error plots
        uErrCumulative{PASS} = zeros(gridHeight, gridWidth);
        vErrCumulative{PASS} = zeros(gridHeight, gridWidth);
        rErrGradCumulative{PASS} = zeros(gridHeight, gridWidth);
        rErrDirectCumulative{PASS} = zeros(gridHeight, gridWidth);
        velocityErrorMagnitudeCumulative{PASS} = zeros(gridHeight, gridWidth);

        % Loop over the sets
        tic

       % ParametersPath
        parametersPath = JobInfo(n).Sets(1).Parameters.Path;

        % Load and parse the image generation job file.
        imageGenerationParameters = load(parametersPath);

        % Number of images generated per set
        imagesPerSet = imageGenerationParameters.JOBFILE.Parameters.Sets.ImagesPerSet;

        % Vortex parameters
        vortexParameters = imageGenerationParameters.VortexParameters;

        % Flow type
        flowType = vortexParameters.VortexType;

        % Flags for lamb vortex or rankine vortex.
        isLamb = ~isempty(regexpi(flowType, 'lam'));
        isRankine = ~isempty(regexpi(flowType, 'ran'));

        % Image times in sort-of physical units (i.e. "seconds")
        imageTimes = imageGenerationParameters.T;

        % Determine the solution times corresponding to the frame numbers
        firstImageTimes  = imageTimes(firstImageNumbers);
        secondImageTimes = imageTimes(secondImageNumbers);

        % Set the solution time to halfway between the image times.
        % This is second-order accurate
        solutionTimes = firstImageTimes + (secondImageTimes - firstImageTimes) / 2;

        % gridpoints
        gridPointsX = (firstVectorFile.X{PASS}(:));
        gridPointsY = (firstVectorFile.Y{PASS}(:));
        
        % Height and width of the interrogation regions for this pass
        regionWidth = JobFile.Parameters.Processing(PASS).InterrogationRegion.Width;
        regionHeight = JobFile.Parameters.Processing(PASS).InterrogationRegion.Height;

        % These are the coordinates of the geometric centroids of each PIV window.
        % This corresponds to the physical location in the flow that the PIV
        % correlation is supposed to measure. For odd sized windows, the geometric
        % centroid of the window is at the center pixel. For even-sized windows, it
        % is 0.5 pixels to the right of the pixel located at (height/2) or
        % (width/2). 
        xCenter_01 = gridPointsX + 0.5 * (1 - mod(regionWidth,  2));
        yCenter_01 = gridPointsY + 0.5 * (1 - mod(regionHeight, 2));
        
        % Vector containing all the x and y grid points
        gridPointsVector = cat(1, gridPointsX, gridPointsY);

        % Physical time elapsed between the two images (i.e. in "seconds")
        interFrameTime = secondImageTimes(1) - firstImageTimes(1);
        
         % Calculate the true solution for lamb vortex flow
        if isLamb
            
            % This reads the vortex parameters corresponding to the 
            % synthetic vortex ring images that were generated
            VortexRadius = vortexParameters.VortexRadius;
            Angle = vortexParameters.Angle;
            R = vortexParameters.CoreRadius;
            XC = vortexParameters.XC;
            YC = vortexParameters.YC;
            xo1 = XC + VortexRadius * sin(Angle);
            xo2 = XC - VortexRadius * sin(Angle);
            yo1 = YC + VortexRadius * cos(Angle);
            yo2 = YC - VortexRadius * cos(Angle);
            Uo = interFrameTime * vortexParameters.PeakVelocity;
            A = 1.25643;
%             
            % This calculates the true Lagrangian velocity field.
            if isLagrangian
                % This calculates the positions of ideal tracer particles after
                % the number of frames defined by "correlationStep"
                [x2, y2] = lambOseenVortexRing(xCenter_01, yCenter_01, [0, correlationStep], vortexParameters);

                % These are the horizontal and vertical components of the
                % analytical Lagrangian velocity field
                uTrueVect = x2(2, :) - x2(1, :);
                vTrueVect = y2(2, :) - y2(1, :);
               
                % This calculates the analytical vorticity and shear
                [~, ~, w_analytical, shear_analytical] = analyticalLambVortex(xCenter_01, yCenter_01, xo1, xo2, yo1, yo2, R, Uo, A);
            else
                % This calculates the true Eulerian velocity field.
                [trueVelocities, trueVort] = lambOseenVortexRingVelocityFunction(0, [xCenter_01; yCenter_01], vortexParameters);
                
                % Extract true displacements
                uTrueVect = interFrameTime * trueVelocities(1:length(trueVelocities)/2);
                vTrueVect = interFrameTime * trueVelocities(length(trueVelocities)/2 + 1 : end);
                rTrueVect = interFrameTime * trueVort(:);

                % Velocity components as matrices
                uTrue{PASS} = flipud((reshape(uTrueVect, size(X{PASS}))));
                vTrue{PASS} = -1 * flipud(reshape(vTrueVect, size(X{PASS})));
                rTrue{PASS} = flipud((reshape(rTrueVect, size(X{PASS}))));
                
                
%                 [uTrueVect, vTrueVect, w_analytical, shear_analytical] = analyticalLambVortex(xCenter_01, yCenter_01, xo1, xo2, yo1, yo2, R, Uo, A);
            end % End if isLagrangian

        else
           % If not a lamb vortex, set the solution to zero for now so things don't break for other flow types.
           trueVelocities = zeros(size(gridPointsVector)); 
           trueVort = zeros(size(gridPointsVector));
        end % End if isLamb

      
%         % This is the true vorticity and shear as a vector.
        % Might not be necessary.
%         rTrueVect = w_analytical(:);
%         shearVect = shear_analytical(:);
% 
        % Velocity components as matrices
        uTrue{PASS} = flipud((reshape(uTrueVect, [gridHeight, gridWidth])));
        vTrue{PASS} = -1 * flipud(reshape(vTrueVect, [gridHeight, gridWidth]));
        rTrue{PASS} = flipud((reshape(rTrueVect, [gridHeight, gridWidth])));
%         shearTrue{PASS} = flipud((reshape(shearVect, [gridHeight, gridWidth])));
        
        
        % Calculate which grid points are inside the vortex "core"
        xInside{PASS} = X{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
        yInside{PASS} = Y{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
        uTrueInside{PASS} = uTrue{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
        vTrueInside{PASS} = vTrue{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
%         trueShearInside{PASS} = shearVect(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
        
        % Number of vectors inside the core
        numInside{PASS} = numel(xInside{PASS});
%         nPoints{PASS} = gridHeight * gridWidth;
        
        % Initialize the file number
        fileNum = 1;
        
        % Number of files to be loaded
        numFiles = nSets * nPairs;
        
        velocityErrorMagnitudeInside{PASS} = zeros(numInside{PASS}, numFiles);
        rErrorGradInside{PASS} = zeros(numInside{PASS}, numFiles);
        rErrorDirInside{PASS} = zeros(numInside{PASS}, numFiles);
        rDirInside{PASS} = zeros(numInside{PASS}, numFiles);
        rGradInside{PASS} = zeros(numInside{PASS}, numFiles);
        spatialPeakRatioInside{PASS} = zeros(numInside{PASS}, numFiles);
        
        % Initialize the vectors to hold all of the errors, from
        % every vector in every pass in ever set. These are big.
        allErr{PASS} = zeros(numFiles * gridHeight * gridWidth, 1);
        allErrInside{PASS} = zeros(numFiles * numInside{PASS}, 1);
        
        % Initialize the vectors to hold the calculated uncertainties
        % calculated from the disparity code.
        displacement_uncertainty_x{PASS} = zeros(numFiles * gridHeight * gridWidth, 1);
        displacement_uncertainty_y{PASS} = zeros(numFiles * gridHeight * gridWidth, 1);
        
        % Initialize vectors to hold the calculated uncertainties near the cores.
        displacement_uncertainty_inside_x{PASS} = zeros(numFiles * numInside{PASS}, 1);
        displacement_uncertainty_inside_y{PASS} = zeros(numFiles * numInside{PASS}, 1);
        
        % Spatial peak ratio vector initialization.
        allSpatialPeakRatio{PASS} = zeros(numFiles * gridHeight * gridWidth, 1);
        allSpatialPeakRatioInside{PASS} = zeros(numFiles * numInside{PASS}, 1);
        
        
    end % End for PASS = 1 : nPasses
    
    % Counter;
    c = 0;

    
    for s = 1 : nSets
%         try
    
        % Inform the user
        disp(['Plotting set ' num2str(s) ' of ' num2str(nSets)])

        % Vector directory
        vectorDir = JobInfo(n).Sets(s).Vector.Directory;
    
        % error directory
        errorDir = JobInfo(n).Sets(s).Error.Directory;
        if ~exist(errorDir, 'dir');
            mkdir(errorDir);
        end

        % Loop over all the image pairs.
        for p = 1 : nPairs
            c = c + 1;
%             disp(['Pair ' num2str(p) ' of ' num2str(nPairs)]);

            % Name of the PIV vector file.
            vectorFileName = [vectorFileBaseName ...
                num2str(firstImageNumbers(p), numberFormat) ...
                '_' num2str(secondImageNumbers(p), numberFormat) '.mat'];
            
            % Path to the PIV vector file
            vectorFilePath = fullfile(vectorDir, vectorFileName);
            
            % Load the PIV output file
            if exist(vectorFilePath, 'file')
                
                % Load the vector file if it exists
                 vectorFile = load(vectorFilePath); 
%                 
%                 if isDeform
%                     % Add to the number of deforms
%                     nDeforms(c, :) = vectorFile.NDEFORMS;
%                 end
                
                % Calculate errors for each pass. 
                for PASS = 1 : nPasses 
                    
                    % Read the PIV-calculated velocities
                    U = vectorFile.U{PASS};
                    V = vectorFile.V{PASS};
                    R = vectorFile.R{PASS};
                    UVAL = vectorFile.UVAL{PASS};
                    VVAL = vectorFile.VVAL{PASS};
                    RVAL = vectorFile.RVAL{PASS};
                    
                    if isfield(vectorFile, 'DISPARITY_X')
                        UNCERT_X = vectorFile.DISPARITY_X{PASS};
                        UNCERT_Y = vectorFile.DISPARITY_Y{PASS};
                    else
                        UNCERT_X = zeros(size(U));
                        UNCERT_Y = zeros(size(U));
                    end
                    
                    
%                     spatialPeakRatio = vectorFile.SPATIAL_PEAK_RATIO{PASS};
            
                    % Number of grid points
                    nPoints = length(U(:));
                    
                    % Number of points inside the cores
                    nPointsInside = numInside{PASS};

                    % Calculate vorticity from direct measurement of rotation.
                    % Multiply 2 because vorticity is 2 * rotation(derivation in notebook)
                    if isLeftHanded
                        rDirect = -(2 * R);
                    else
                        rDirect =  (2 * R);
                    end

                    % Calculate vorticity from gradient of velocity field.
                    % Calculate gradients using second-order central differencing scheme
                    % and use the validated vector field. Don't validate
                    % the vectors in-line because that takes a while;
                    % instead, this should have been done during vector
                    % calculation.
                    rGrad = calculateVorticity_socdiff(X{PASS}, Y{PASS}, UVAL, VVAL, RVAL, 0);
                                    
                    % Calculate velocity errors (signed)
                    uErr = (U - uTrue{PASS});
                    vErr = (V - vTrue{PASS});

                    % Calculate velocity error magnitude
                    absVelocityError = sqrt(uErr.^2 + vErr.^2);
                   
                   % Error of direct vorticity estimate
                    rErrDirect  = abs(rTrue{PASS} - rDirect);
                    
                   % Error of gradient-based vorticity estimate
                    rErrGrad = abs(rTrue{PASS} - rGrad);
               
                    % Calculate the error metrics near the core
                    errInside = absVelocityError(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    rErrGradInside_pass = rErrGrad(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    rErrDirInside_pass = rErrDirect(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    uncert_x_inside = UNCERT_X(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    uncert_y_inside = UNCERT_Y(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    meanErrInside(c, PASS) = mean(errInside(:));
                    
%                     spatialPeakRatioInside = spatialPeakRatio(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    
                    firstPointAll = sub2ind([nPoints, nPairs, nSets], 1, p, s);
                    lastPointAll = sub2ind([nPoints, nPairs, nSets], nPoints, p, s);
                    
                    % Determine points inside the core
                    firstPointInside = sub2ind([nPointsInside, nPairs, nSets], 1, p, s);
                    lastPointInside =  sub2ind([nPointsInside, nPairs, nSets], nPointsInside, p, s);
                    
                    % Velocity error magnitudes
                    allErr{PASS}(firstPointAll:lastPointAll) = absVelocityError(:);
                    allErrInside{PASS}(firstPointInside:lastPointInside) =  errInside(:);
                    
                    % Displacement uncertainties 
                    displacement_uncertainty_x{PASS}(firstPointAll:lastPointAll) = UNCERT_X(:);
                    displacement_uncertainty_y{PASS}(firstPointAll:lastPointAll) = UNCERT_Y(:);
                    
                    % Displacement uncertaintites near the core
                    displacement_uncertainty_inside_x{PASS}(firstPointInside:lastPointInside) = uncert_x_inside(:);
                    displacement_uncertainty_inside_y{PASS}(firstPointInside:lastPointInside) = uncert_y_inside(:);
                    
%                     allSpatialPeakRatio{PASS}(firstPointAll:lastPointAll) = spatialPeakRatio(:);
%                     allSpatialPeakRatioInside{PASS} = spatialPeakRatioInside(:);

                    
                    
                    

                    velocityErrorMagnitudeInside{PASS}(:, fileNum) = errInside(:);
                    rErrorGradInside{PASS}(:, fileNum) = rErrGradInside_pass(:);
                    rErrorDirInside{PASS}(:, fileNum) = rErrDirInside_pass(:);
%                     spatialPeakRatioInside{PASS}(:, fileNum) = peakRatioInside(:);

                    
                    rdi = rDirect(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    rgi = rGrad(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
                    rDirInside{PASS}(:, fileNum) = rdi(:);
                    rGradInside{PASS}(:, fileNum) = rgi(:);
                    
                    


    %                 imagesc(X(:), Y(:), velocityErrorMagnitude);
    %                 caxis([0 3]);
    %                 hold on;
    %                 quiver(X, Y, uTrue, vTrue, 'white');
    %                 quiver(X, Y, U, V, 'black');
    %                 plot(xInside, yInside, '.w');
    %                 axis image;
    %                 hold off
    %                 set(gca, 'YDir', 'Reverse');
    %                 title(['Set ' num2str(s) ', Pair ' num2str(p)], 'FontSize', 18);
    %                 drawnow

                    % Cumulative errors
                    uErrCumulative{PASS} = uErrCumulative{PASS} + abs(uErr);
                    vErrCumulative{PASS} = vErrCumulative{PASS} + abs(vErr);
                    rErrGradCumulative{PASS} = rErrGradCumulative{PASS} + abs(rErrGrad);
                    rErrDirectCumulative{PASS} = rErrDirectCumulative{PASS} + abs(rErrDirect);
                    velocityErrorMagnitudeCumulative{PASS} = velocityErrorMagnitudeCumulative{PASS} + absVelocityError;

                end % end for PASS = 1 : nPasses
                fileNum = fileNum + 1;
            end % end if exist(vectorFilePath, 'file');
  
        end % end for p = 1 : nPairs
%         catch
% %             keyboard
%         end
    end
    
    % Average number of deform passes 
    %nDeformsAverage = mean(nDeforms, 1);
   % nDeformsStd = std(nDeforms, 1);
    
    
    toc;

    % Calculate the average errors for each pass.
    for PASS = 1 : nPasses
        % Average the errors.
        uErrAve{PASS} = uErrCumulative{PASS} ./ (nPairs * nSets);
        vErrAve{PASS} = vErrCumulative{PASS} ./ (nPairs * nSets);
        rErrGradAve{PASS} = rErrGradCumulative{PASS} ./ (nPairs * nSets);
        rErrDirectAve{PASS} = rErrDirectCumulative{PASS} ./ (nPairs * nSets);
        velocityErrorMagnitudeAve{PASS} = velocityErrorMagnitudeCumulative{PASS} ./ (nPairs * nSets);
        
% Calculate the average errors inside the vortex core.
%         uErrInside{PASS} = uErr(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:))); 
%         vErrInside{PASS} = vErr(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
%         velocityErrorMagnitudeInside{PASS} = velocityErrorMagnitudeAve{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
%         rErrGradInside{PASS} = rErrGradAve{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
%         rErrDirInside{PASS} = rErrDirectAve{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
        rTrueInside{PASS} = rTrue{PASS}(abs(rTrue{PASS}) >= vorticityThreshold * max(rTrue{PASS}(:)));
    end
    
    % Extract some other information about the PIV job
     imageBaseName = JobFile.Parameters.Images.BaseName;
     particleConcentration = JobFile.Parameters.Images.ParticleConcentration;
     startSet = JobInfo(n).Sets(1).SetNumber;
     endSet   = JobInfo(n).Sets(end).SetNumber;
     
     % Specify name of the error statistics file
     errorName = ['errorStats_' correlationMethod '_'...
         imageBaseName  'c_' num2str(particleConcentration, '%0.4f')...
         '_cstep_' num2str(correlationStep, '%02.0f') '_sets_' num2str(startSet, '%05.0f') '-' num2str(endSet, '%05.0f') '.mat'];
     
     % Specify path to the error statistics file
     errorPath = fullfile(statsDir, errorName);
     
     % Save the error file
        save(errorPath, 'X', 'Y', ...
            'uTrue', 'uTrueInside', 'vTrueInside',...
            'vTrue', 'rTrue', 'rTrueInside', 'rGradInside',...
            'rDirInside', 'xInside', 'yInside',...
            'uErrAve', 'vErrAve', 'velocityErrorMagnitudeAve',...
            'velocityErrorMagnitudeInside', 'rErrDirectAve',...
            'rErrorDirInside', 'rErrGradAve', 'rErrorGradInside',...
            'vorticityThreshold', 'allErr', 'allErrInside',...
            'allSpatialPeakRatio', 'allSpatialPeakRatioInside', ...
            'displacement_uncertainty_x', 'displacement_uncertainty_y', ...
            'displacement_uncertainty_inside_x', 'displacement_uncertainty_inside_y', 'meanErrInside');
end

    


end


