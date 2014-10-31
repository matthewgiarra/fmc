function IMAGE = im2logpol(IMAGE, NUMWEDGES, NUMRINGS)
% im2logpol(IMAGE) calculates the log-polar transform of an image and
% outputs the result as an 8-bit image

% % Load the image
% Read the image if the input is a string referring to a file path
if ischar(IMAGE)
    img = imread(IMAGE);
else
    img = IMAGE;
end

% Determine size of image. nChannels is the number of color channels
[height, width, nChannels] = size(img);

% Coordinates of center of image
xCenter = (width + 1) / 2; 
yCenter = (height + 1) / 2;

% Specify min and max radius of image
% maxRadius = sqrt(xCenter^2 + yCenter^2);
% maxRadius = floor(min(width/2, height/2));

% % Number of circular sections
% if nargin < 3
%     nRings = maxRadius;
% %         nRings = width;
% else
%     nRings = NUMRINGS;
% end
% 
% % Number of angular sections
% if nargin < 2
%     nWedges = ceil(2 * pi * maxRadius);
% %        nWedges = height;
% else
%     nWedges = NUMWEDGES;
% end
% 
% if mod(nWedges, 2) > 0
%     nWedgesEven = nWedges + 1;
% else
%     nWedgesEven = nWedges;
% end


% % TESTING!!
% nWedges = height;
% nRings = width;

% nWedges = 2000;
% nWedges = maxRadius;

% Spatial location of the image in the 2-D input space
UData = [1, width] - xCenter;
VData = [1, height] - yCenter;

% % Spatial location of the image in the 2-D output space
% XData = [0, log(maxRadius)];
% YData = [0, nWedgesEven - 1];

% Spatial location of the image in the 2-D output space
XData = [0, log(NUMRINGS)];
YData = [0, NUMWEDGES - 1];

% Number of rows and columns of the output image
% Size = [height, width];
Size = [NUMWEDGES, NUMRINGS];

% Create image transformation structure
t = logpoltform(NUMWEDGES);

% Transform a single-channel image
if nChannels == 1
    IMAGE = imtransform(img, t, 'bicubic', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'Size', Size); 
else
        IMAGE = zeros(height, width, nChannels); % Initialize the output image
        for k = 1:nChannels % Loop through all the channels
                IMAGE(:, :, k) = imtransform(img(:, :, k), t, 'bicubic', 'UData', UData, 'VData', VData, 'XData', XData, 'YData', YData, 'Size', Size); % Transform each channel
        end        
end

% Set infinities to zero and save information to the output variable. This
% is necessary when the input image is in the frequency domain.
IMAGE(isinf(IMAGE)) = 0;

end
