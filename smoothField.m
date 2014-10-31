function [VECTOR_OUT] = smoothField(VECTOR_IN, GAUSSIAN_SIZE, GAUSSIAN_STD)

% Create the gaussian filter.
gaussianFilter = fspecial('gaussian', GAUSSIAN_SIZE, GAUSSIAN_STD);

%there seems to be a very strange bug in imfilter for when U is single and
%small in size (<18?).  Casts are very fast, and forcing it to double fixes 
%the problem.
VECTOR_OUT = cast(imfilter(double(VECTOR_IN), gaussianFilter,'replicate'), class(VECTOR_IN));


end