function FILEPATHS = makeFilepaths(IMAGEDIRECTORY, IMAGEBASENAME, NUMBERFORMAT, FILEEXTENSION, STARTFILE, STOPFILE, FILETRAILER, IMAGESTEP)
%  MAKEILEPATHS creates a vector whose rows contain the file-paths to
%  sequentially-numbered files
% 

% % Read variables from structure
% imDir = INPUTS.Inputs.ImDir;  % Directory containing raw images
% imBase = INPUTS.Inputs.ImBase; % Basename of raw images
% nDigits = INPUTS.Inputs.nDigits; % Number of digits in file numbers
% fileExtension = INPUTS.Inputs.FileExtension; % Extension
% startFile = INPUTS.Inputs.StartImage; % First image
% stopFile= INPUTS.Inputs.StopImage;  % Last image
% fileTrailer = INPUTS.Inputs.FileTrailer; % Text trailing the file number but preceeding the extension
% imageStep = INPUTS.Options.ImageStep; % Step between images

depth = STOPFILE - STARTFILE + 1; % Number of images

% Create array of strings containing paths to
% for k = 1 : IMAGESTEP : depth
%    filepaths(k, :) = fullfile(IMAGEDIRECTORY, [ IMAGEBASENAME sprintf(NUMBERFORMAT, STARTFILE + k - 1) FILETRAILER FILEEXTENSION]);
% end

for k = 1 : 1 : floor(depth/IMAGESTEP)
   filepaths(k, :) = fullfile(IMAGEDIRECTORY, [ IMAGEBASENAME sprintf(NUMBERFORMAT, STARTFILE + (k - 1) * IMAGESTEP) FILETRAILER FILEEXTENSION]);
end


FILEPATHS = filepaths; % Save to output variable

end