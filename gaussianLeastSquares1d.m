function GAUSSIANCENTROIDLOCATIONOUTPUT = gaussianLeastSquares1d(ARRAY, WINDOWRADIUS)

%%%% TEST VARIABLES. DELETE WHEN FUNCTIONALIZED.

% % Load the test data
% Data = load('~/Desktop/spatialCorr.mat');
% ARRAY = Data.spatialCorr;

% Radius of the window over which to fit the Gaussian function.
% The input data used in the fitting will be symmetric about the
% maximum value of the array, and the total size of the array used 
% in the fit will be 2 * WINDOWRADIUS + 1
% WINDOWRADIUS = 2;

%%%% END TEST VARIABLES.

% Length of the signal (number of elements in the array)
arrayLength = length(ARRAY);

% Find the magnitude and the location of the maximum value of the array.
% This should be updated to deal with multiple values that are equal to one another. 
% maxValue is the maximum value in the signal, and maxIndex is the index
% location of the max value within the data.
[arrayMax, arrayMaxIndex] = max(ARRAY);

% Determine the minimum and maxium indices of the data to which
% the gaussian function be fit. If the window exceeds the extent
% of the array, terminate the window at the beginning or the end of the
% array.
minIndex = max(1, arrayMaxIndex - WINDOWRADIUS);
maxIndex = min(arrayLength, arrayMaxIndex + WINDOWRADIUS);

% Form a vector containing the indices of the data to which
% the Gaussian function will be fit.
fitArrayIndices = (minIndex : maxIndex)';

% Form a vector containing the values of the data to which
% the Gaussian function will be fit.
fitArrayData = ARRAY(fitArrayIndices);

% Specify initial guesses for the parameters describing the Gaussian
% function that will be fit to the data. The parameters that will be solve
% for are (1) the max value of the Gaussian function (gaussianMaxValue), (2) the location of
% the centroid (gaussianCentroidLocation) of the Gaussian function, and (3) the standard deviation of
% the Gaussian function (gaussianStdDev). 

% Guess that the max value of the Gaussian function is
% equal to the max value of the input array.
gaussianMaxValueGuess = arrayMax;

% Guess that the location of the centroid of the Gaussian function
% is equal to the location of the max value of the input array.
gaussianCentroidLocationGuess = arrayMaxIndex;

% Guess that the standard deviation of the Gaussian function
% is equal to one.
gaussianStdDevGuess = 1;

% Guess the noise floor of the Gaussian function is equal to zero
gaussianNoiseFloorGuess = 0;

% Make an array containing the initial guesses for the parameters of
% the Gaussian function. This array will be passed to the least squares solver.
gaussianParametersGuess = [gaussianMaxValueGuess, ...
    gaussianCentroidLocationGuess, gaussianStdDevGuess, gaussianNoiseFloorGuess];

% Set options for the least squares solver.
solverOptions = optimset('MaxIter', 1200, 'MaxFunEvals', 5000,...
    'TolX', 5e-6,'TolFun', 5e-6, 'LargeScale', 'off', 'Display', 'off', ...
    'DiffMinChange', 1e-7, 'DiffMaxChange', 1, 'Algorithm','levenberg-marquardt');

% Set the lower and upper bounds of the solutions to be empty arrays
% so that the bounded solver is not used. The purpose of coding the
% function this way is so that it can be modified to calculate
% bounded solutions if desired.
solutionLowerBounds = [];
solutionUpperBounds = [];

% Run the least squares solver and return a vector that contains the
% parameters of the best-fit Gaussian function. The elements
% gaussianParametersOutput are:
% gaussianParametersOutput(1) = max value of the best-fit Gaussian function
% gaussianParametersOutput(2) = location of the centroid of the best-fit Gaussian function
% gaussianParametersOutput(3) = standard deviation of the best-fit Gaussian function
% gaussianParametersOutput = lsqnonlin(@leastsquares1D, gaussianParametersGuess, ...
%     solutionLowerBounds, solutionUpperBounds, solverOptions, fitArrayData, fitArrayIndices);
gaussianParametersOutput = lsqnonlin(@gaussianDifference1D, gaussianParametersGuess, ...
    solutionLowerBounds, solutionUpperBounds, solverOptions, fitArrayIndices, fitArrayData );

% Save the location of the centroid of the best-fit Gaussian
% to the output variable.
GAUSSIANCENTROIDLOCATIONOUTPUT = gaussianParametersOutput(2);

% % Generate a plot comparing the input data to the best fit Gaussian
% % Uncomment for debugging.
% gaussianMaxValueOutput = gaussianParametersOutput(1);
% gaussianVarianceOutput = gaussianParametersOutput(3);
% gaussianNoiseFloorOutput = gaussianParametersOutput(4);
% x = minIndex : 0.01 : maxIndex;
% y = gaussianMaxValueOutput * ...
%     exp(-(x - GAUSSIANCENTROIDLOCATIONOUTPUT).^2 / gaussianVarianceOutput ) ...
%     + gaussianNoiseFloorOutput;
% % 
% plot(ARRAY);
% hold on
% plot(fitArrayIndices, fitArrayData, 'or');
% plot(x, y, '--k');
% hold off
% xlabel('Index', 'FontSize', 16);
% ylabel('Value', 'FontSize', 16);
% title('Comparison of input data and best-fit Gaussian function', 'FontSize', 16);
% set(gca, 'FontSize', 16);
% set(gcf, 'color', 'white');
% legend('Input data', 'Best-fit Gaussian', 'FontSize', 16);

end







