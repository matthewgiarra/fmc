function [XLP, YLP] = LogPolarCoordinates(IMAGESIZE, NUMWEDGES, NUMRINGS, RMIN, RMAX, MAXANGLE)

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
% xZero = w/2 + 1 - 0.5 * mod(w, 2);
% yZero = h/2 + 1 - 0.5 * mod(h, 2);

xZero = (w + 1)/2;
yZero = (h + 1)/2;

% Log(r) coordinate vector
logR = linspace(log(rMin), log(rMax), nr);

% Radial coordinate vector
rv = exp(logR);

% Angular coordinate vector
thMax =  maxAngle * (1 - 1 / nw);
% thMax =  maxAngle;

thv = linspace(0, thMax, nw);

% Log-polar grid
[r, th] = meshgrid(rv, thv);

% Convert log polar grid to cartesian.
[x, y] = pol2cart(th, r);


% [r1, th1] = meshgrid(rv, thv(1:nw/2 + 1));
% [x1, y1] = pol2cart(th1, r1);
% 
% [r2, th2] = meshgrid(rv, thv(nw/2+2 : end));
% [x2, y2] = pol2cart(th2, r2);

% Shift the center of the log-polar grid to the zero-frequency
% element of the image. Assign these coordinate matrices
% to the output variables.
XLP = x + xZero;
YLP = y + yZero;

% XLP1 = x1 + 1;
% YLP1 = y1 + 1;
% 
% XLP2 = x2 + 2 * (xZero - 1);
% YLP2 = y2 + 1;
% 
% XLP = cat(1, XLP1, XLP2);
% YLP = cat(1, YLP1, YLP2);

end
