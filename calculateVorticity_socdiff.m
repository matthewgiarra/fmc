% function [CIRCULATIONGRADIENTMETHOD CIRCULATIONDIRECTMETHOD LCI UVAL VVAL RVAL VORTICITYGRAD LCIVAL] = ...
%     calculateVorticity(X, Y, U, V, R, THRESHOLD, VALIDATE)

function [VORTICITYGRAD, UVAL, VVAL, RVAL] = calculateVorticity_socdiff(X, Y, U, V, R, VALIDATE)

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
    UODwinsize = [3 3; 3 3];
    BootPer = 0;
    BootIter = 0;
    BootKMax = 0;

    % Vector validation
    [UvalRaw,VvalRaw, RvalRaw] = VAL_2d3c(X(:),Y(:),U(:),V(:), R(:), Eval(:),C(:),D(:),ThreshSwitch,UODswitch,...
        BootSwitch, extraPeakSwitch,Uthresh,Vthresh,UODwinsize,UODthresh,BootPer,BootIter,BootKMax);
    UVAL = flipud(UvalRaw);
    VVAL = flipud(VvalRaw);
    RVAL = flipud(RvalRaw);
    
else
     UVAL = U;
     VVAL = V;
     RVAL = R;
end

% Spacing between grid points
dx = X(1,2) - X(1,1);
dy = Y(2, 1) - Y(1, 1);

% Partial derivatives of velocity
dvdx = socdiff(V, dx, 2);
dudy = socdiff(U, dy, 1);

VORTICITYGRAD = -1 * (dvdx - dudy);


end

