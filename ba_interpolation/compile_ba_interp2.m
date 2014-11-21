function compile_ba_interp2
% This function compiles the c++ function ba_interp2.cpp
% into a Matlab mex file.

% Compile ba_interp2
mex -O ../ba_interpolation/ba_interp2.cpp

end