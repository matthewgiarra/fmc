function [SIMILARITY_MATRIX_01, SIMILARITY_MATRIX_02] = makeSimilarityMatrices_homogeneous(ROTATIONANGLE, SCALING, TRANSLATION_Y, TRANSLATION_X, DIFFERENCEMETHOD)
% [SIMILARITY_MATRIX_01, SIMILARITY_MATRIX_02] = makeSimilarityMatrices(ROTATIONANGLE, SCALING, DIFFERENCEMETHOD)
% This function creates a pair of similarity matrices for aligning a pair
% of images in either the central, forward, or backward difference sense.
%
% INPUTS


% Rename the rotation angle for readability
th = ROTATIONANGLE; 

% Rename the input scaling variable for readability
s = SCALING;

% Rename the translations for readability
ty = TRANSLATION_Y;
tx = TRANSLATION_X;

% Differencing methods:
% DIFFERENCEMETHOD = 1 is central difference
% DIFFERENCEMETHOD = 2 is forward difference
% DIFFERENCEMETHOD = 3 is backward difference

switch DIFFERENCEMETHOD
    
    % This is the case of central difference method FMC
    case 1       
        % Rotation matrices 
        % This matrix rotates the FIRST image by half the originally-calculated 
        % rotation angle in the FORWARD-difference direction.
        R1 = [cos( th/2 ) -sin( th/2 ) 0; sin( th/2 ) cos( th/2 ) 0; 0 0 1]; % Matrix originally calculated rotation angle

        % This matrix rotates the SECOND image by half the originally-calculated
        % rotation angle in the BACKWARD-difference direction.
        R2 = [cos( -th/2 ) -sin( -th/2 ) 0; sin( -th/2 ) cos( -th/2 ) 0; 0 0 1];

        % This matrix scales the first image by half (or really square root of) the
        % originally-calculated scaling factor in the FORWARD difference direction.
        S1 = sqrt([ s 0 0; 0 s 0; 0 0 1]);

        % This matrix scales the second image by the inverse of half
        % (or really square root of) the originally-calculated scaling factor
        % in the BACKWARD difference direction.
        S2 = sqrt([1/s 0 0; 0 1/s 0; 0 0 1]);
        
        % This matrix translates the first image by half the calculated
        % translation in the forward difference direction
        T1 = [1 0 tx/2; 0 1 ty/2; 0 0 1];
        
        % This matrix translates the second image by half the calculated
        % translation in the backward difference direction
        T2 = [1 0 -tx/2; 0 1 -ty/2; 0 0 1];
        
    % This is the case of forward-difference FMC
    case 2
        % Rotation matrices 
        % This matrix rotates the FIRST image by the originally-calculated 
        % rotation angle in the FORWARD-difference direction.
        R1 = [cos( th ) -sin( th ) 0; sin( th ) cos( th ) 0; 0 0 1];

        % This leaves the second image untouched.
        R2 = [1 0 0; 0 1 0; 0 0 1];

        % This matrix scales the first image by the originally-calculated
        % scaling factor in the FORWARD difference direction.
        S1 = [ s 0 0; 0 s 0; 0 0 1];

        % This matrix leaves the second image untouched.
        S2 = [1 0 0; 0 1 0; 0 0 1];
        
        % This matrix translates the first image by the originally calculated
        % translations
        T1 = [1 0 tx; 0 1 ty; 0 0 1];
        
        % This matrix leaves the second image untouched.
        T2 = [1 0 0; 0 1 0; 0 0 1];
       
    % This is the case of backward-difference FMC
    case 3        
        % Rotation matrices   
        % This leaves the first image untouched.
        R1 = [1 0 0; 0 1 0; 0 0 1];
        
        % This matrix rotates the SECOND image by the negative of the
        % originally-calculated rotation angle in the BACKWARD-difference direction.
        R2 = [cos( -th ) -sin( -th ) 0; sin( -th ) cos(-th) 0; 0 0 1];
        
        % This leaves the first image untouched.
        S1 = [1 0 0; 0 1 0; 0 0 1];

        % This matrix scales the second image by the inverse of the 
        % originally-calculated scaling factor in the BACKWARD difference direction.
        S2 = [ 1/s 0 0; 0 1/s 0; 0 0 1];
        
        % This matrix leaves the first image untouched
        T1 = [1 0 0; 0 1 0; 0 0 1];
        
        % This matrix translates the second image by the negative of the
        % originally calculated translations.
        T2 = [1 0 -tx; 0 1 -ty; 0 0 1];        
end

% Multiply the translation, scaling, and rotation matrices 
% to get the homogeneous similarity matrices.
SIMILARITY_MATRIX_01 = T1 * S1 * R1;
SIMILARITY_MATRIX_02 = T2 * S2 * R2;

end

