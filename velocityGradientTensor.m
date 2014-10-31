function [VELOCITYGRADIENTTENSOR, SYMMETRICTENSOR, ANTISYMMETRICTENSOR] = velocityGradientTensor(X, Y, U, V)

% Calculate the grid spacing in the X (column) direction
dx =unique(diff(X, 1, 2));
DX = dx(1);

% Calculate the grid spacing in the Y (row) direction
dy = unique(diff(Y, 1, 1)); 
DY = dy(1);

% Calculate the spatial gradients of velocity using the Richardson 4th
% order gradient scheme
[dudx, dudy, dvdx, dvdy] = richardsonGradient(U, V, DX, DY);

% Permute the different components of the spatial velocity gradients so
% that the first two dimensions correspond to the rows and columns of the
% gradient tensor, and the second two dimensions correspond to the rows and
% columns of the vector field. 
dUdX = permute(dudx, [ 3, 4, 1, 2 ] );
dUdY = permute(dudy, [ 3, 4, 1, 2 ] );
dVdX = permute(dvdx, [ 3, 4, 1, 2 ] );
dVdY = permute(dvdy, [ 3, 4, 1, 2 ] );

% Build the velocity gradient tensor
VELOCITYGRADIENTTENSOR = [dUdX dUdY; dVdX dVdY];

% Calculate the symmetric part of the velocity gradient tensor. The
% operation "permute( X, [ 2, 1, 3, 4] )" switches the first two dimensions
% of X, so that " X + permute(X, [ 2, 1, 3, 4] ) = X_i, j + X_j, i ". 
SYMMETRICTENSOR = 1 / 2 * ( VELOCITYGRADIENTTENSOR + permute(VELOCITYGRADIENTTENSOR, [ 2, 1, 3, 4 ]  ) ); 

% Calculate the antisymmetric part of the velocity gradient tensor
ANTISYMMETRICTENSOR = 1 / 2 * ( VELOCITYGRADIENTTENSOR - permute(VELOCITYGRADIENTTENSOR, [ 2, 1, 3, 4 ] ) );


end