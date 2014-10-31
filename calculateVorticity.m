% function [CIRCULATIONGRADIENTMETHOD CIRCULATIONDIRECTMETHOD LCI UVAL VVAL RVAL VORTICITYGRAD LCIVAL] = ...
%     calculateVorticity(X, Y, U, V, R, THRESHOLD, VALIDATE)

function [UVAL VVAL RVAL VORTICITYGRAD LCIVAL] = calculateVorticity(X, Y, U, V, R, VALIDATE)


% THRESHOLD = 0.07;
% THRESHOLD = 0.00;

% Calculate the size of the data (number of rows, number of columns).
[nRows nColumns] = size(X);

if VALIDATE
    Eval = zeros(size(U));
    C = zeros(size(U));
    D = zeros(size(U));
    ThreshSwitch = 0;
    UODswitch = 1;
    BootSwitch = 0;
    extraPeakSwitch = 0;
    Uthresh = 32;
    Vthresh = 32;
    UODthresh = [3 2];
    UODwinsize = [ 3 3; 3 3];
    BootPer = 0;
    BootIter = 0;
    BootKMax = 0;

    % Vector validation
    [UvalRaw,VvalRaw, RvalRaw ~, ~, ~] = VAL_rotation(X(:),Y(:),U(:),V(:),R(:), Eval(:),C(:),D(:),ThreshSwitch,UODswitch,...
        BootSwitch, extraPeakSwitch,Uthresh,Vthresh,UODwinsize,UODthresh,BootPer,BootIter,BootKMax);
    UVAL = flipud(UvalRaw);
    VVAL = flipud(VvalRaw);
    RVAL = flipud(RvalRaw);
    
    
else
     UVAL = U;
     VVAL = V;
     RVAL = R;
end


% Calculate the magnitudes of the imaginary parts of the complex
% eigenvectors of the velocity gradient tensor (Lambda-ci vortex ID method). 
% Points where Lambda-ci > 0 are taken to be inside a vortex. Ref: J. Zhou,
% "Mechanisms for Generating Coherent Packets of Hairpin Vortices in
% Channel Flow." Journal of Fluid Mechanics, 1999
[LCI, velocityGradient, ~, ~] = lambdaci(X, Y, UVAL, VVAL);

% [LciRaw, velocityGradientRaw, ~, ~] = lambdaci(X, Y, U, V);

% Lci squared
Lci2 = LCI.^2;

% Max of Lci squared
maxLci2 = max(Lci2(:));

% Rename the velocity gradient tensor as G for brevity in coding.
G = velocityGradient;

% Find the row and column locations of grid points that are inside a vortex
% (i.e. grid points where Lci > 0)
% [vortexRows, vortexCols] = find(Lci2 >= THRESHOLD * maxLci2);
% LCIVAL = Lci2 >= THRESHOLD * maxLci2;

% Convert the row and column positions to a single-digit linear index.
% vortexIndices = sub2ind(size(LCI), vortexRows, vortexCols);

% Initialize the vorticity matrix
VORTICITYGRAD = zeros(nRows, nColumns);

% Calculate the vorticity in the flow field
for m = 1 : nRows
    for n = 1 : nColumns
      VORTICITYGRAD(m, n) =   - 1 / 2 * ( G(2, 1, m, n) - G(1, 2, m, n));
    end
end

% % Calculate circulation by adding up the vorticity over all of the points
% % inside of the vortices (i.e. surface integral of vorticity within the vortices). 
% % CIRCULATION = sum(vorticity(vortexIndices));
% % if strcmpi(METHOD, 'direct')
%     CIRCULATIONDIRECTMETHOD = sum(abs(RVAL(vortexIndices)));
% % else
%     CIRCULATIONGRADIENTMETHOD = sum(abs(VORTICITYGRAD(vortexIndices)));
% % end



% contourf(X, Y, Lci2 >= THRESHOLD * maxLci2);
% axis image
% hold on;
% quiver(X, Y, U, V, 'w');
% hold off
% pause(0.01);


end

