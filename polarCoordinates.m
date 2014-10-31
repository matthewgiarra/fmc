function [XLP, YLP] = polarCoordinates(IMAGESIZE, NUMWEDGES, NUMRINGS, RMIN, RMAX, MAXANGLE)

if nargin < 6   
    maxAngle = pi;
else
    maxAngle = MAXANGLE;
end

% Grid dimensions
h = IMAGESIZE(1);
w = IMAGESIZE(2);

nw = NUMWEDGES;
nr = NUMRINGS;
rMax = RMAX;
rMin = RMIN;

% Calculate the index in the array that corresponds
% to the zero frequency. For dimensions with an odd
% number of elements, the zero frequency along that 
% dimension is located at exactly the center of the 
% dimension, i.e., at length / 2 + 0.5. For dimensions
% with an even number of elements, the zero frequency
% is located one pixel to the right of the center, i.e.,
% at length / 2 + 1. 
xZero = w/2 + 1 - 0.5 * mod(w, 2);
yZero = h/2 + 1 - 0.5 * mod(h, 2);

% Log(r) coordinate vector
% logR = linspace(log(rMin), log(rMax), nr);
% Radial coordinate vector
% rv = exp(logR);

rv = linspace(rMin, rMax, nr);

% Angular coordinate vector
thMax =  maxAngle * (1 - 1 / nw);
thv = linspace(0, thMax, nw);

% Log-polar grid
[r, th] = meshgrid(rv, thv);

% Convert log polar grid to cartesian.
[x, y] = pol2cart(th, r);

% Shift the center of the log-polar grid to the zero-frequency
% element of the image. Assign these coordinate matrices
% to the output variables.
XLP = x + xZero;
YLP = y + yZero;



end
