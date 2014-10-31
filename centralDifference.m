function DIFF = centralDifference(MATRIX, DIMENSION)

% Transpose the matrix so that the desired differencing direction is
% column-wise
if DIMENSION == 1
    paddedData = padarray(MATRIX, [1 1], 0);
elseif DIMENSION == 2;
    paddedData = (padarray(MATRIX, [1 1], 0))';
end

% Size of the padded array
[ nRows nCols ] = size(paddedData);

% Number of elements in the padded array
nElements = nRows * nCols;

% Reshape the data matrix into a column vector so that the column is
% aligned with the desired differencing direction.
dataVector = reshape(paddedData, nElements, 1);

% Initialize the difference vector 
diffVect = zeros(nElements, 1);
diffVect(1) = 0;
diffVect(end) = 0;

for k = 2 : nElements - 1
    diffVect( k ) = 1 / 2 * (dataVector( k + 1) - dataVector( k - 1));
end

% Reshape the difference vector back into a matrix
diffMat = reshape(diffVect, nRows, nCols);

% Decide whether or not to re-transpose the matrix
if DIMENSION == 1
    paddedDiff = diffMat;
elseif DIMENSION == 2
    paddedDiff = diffMat';
end

% Strip off the zero-padded part of the matrix
DIFF = paddedDiff(2 : end - 1, 2 : end - 1);

end





