function [UVAL, VVAL, WVAL, ISOUTLIER] = validateField_2d3c(X, Y, U, V, W)

% Hard code a bunch of validation parameters
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
[UvalRaw, VvalRaw, RvalRaw, ISOUTLIER] = VAL_2d3c(X(:), Y(:), U(:), V(:), W(:), Eval(:), C(:), D(:), ThreshSwitch,UODswitch,...
    BootSwitch, extraPeakSwitch, Uthresh, Vthresh, UODwinsize, UODthresh, BootPer,BootIter,BootKMax);
UVAL = flipud(UvalRaw);
VVAL = flipud(VvalRaw);
WVAL = flipud(RvalRaw);
   
end

