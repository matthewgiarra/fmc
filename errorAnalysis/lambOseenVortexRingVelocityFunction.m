
function [U, VORTICITY] = lambOseenVortexRingVelocityFunction(T, X, VORTEXPARAMETERS)

% circulation = VORTEXPARAMETERS.Circulation;
vortexRadius = VORTEXPARAMETERS.VortexRadius;
coreRadius = VORTEXPARAMETERS.CoreRadius;
XC = VORTEXPARAMETERS.XC;
YC = VORTEXPARAMETERS.YC;
vortexAngle = VORTEXPARAMETERS.Angle;
% viscosity = VORTEXPARAMETERS.Viscosity;
propagationVelocity = VORTEXPARAMETERS.PropagationVelocity;
maxTangentialVelocity = VORTEXPARAMETERS.PeakVelocity;

% % Convert input coordinates to column vectors
% x = X(1, :);
% y = X(2, :);

% Number of points
nPoints = numel(X)/2;

% Extract the horizontal and vertical coordinates from the column-vector
% input. The input is a column vector because it has to work with ODE45,
% which demands a column vector input.
x = X(1 : nPoints);
y = X(nPoints + 1 : end);

% Center of the vortex ring propagated by the vortex velocitity
xc = XC + propagationVelocity .* T * cosd(vortexAngle);
yc = YC + propagationVelocity .* T * sind(vortexAngle);

% Positions of vortex cores (first core)
xc1 = xc - vortexRadius * sind(vortexAngle);
yc1 = yc + vortexRadius * cosd(vortexAngle);

% Positions of vortex cores (second core)
xc2 = xc + vortexRadius * sind(vortexAngle);
yc2 = yc - vortexRadius * cosd(vortexAngle);

% Transform the particle coordinates to polar coordinates
[th1, r1] = cart2pol(x - xc1 + 1, y - yc1 + 1);
[th2, r2] = cart2pol(x - xc2 + 1, y - yc2 + 1);

% % Number of positions (i.e. particles)
% nPositions = numel(x);

% % Initialize the angular positions of the particles in the second image
% uTheta1 = zeros(nPoints, 1);
% uTheta2 = zeros(nPoints, 1);

% Tangential velocities
% uTheta1 =       circulation ./ (2 * pi * r1) .* (1 - exp(-r1.^2 ./ (4 * viscosity * T)));
% uTheta2 = - 1 * circulation ./ (2 * pi * r2) .* (1 - exp(-r2.^2 ./ (4 * viscosity * T)));

% Parameter (alpha) defined by Davenport et al, 1996, "The
% structure and development of a wing-tip vortex"
alpha = 1.25643;

uTheta1 =      maxTangentialVelocity * (1 + 0.5 / alpha ) * (coreRadius ./ r1) .* (1 - exp(-1 * alpha  .* r1.^2 / coreRadius^2));
uTheta2 = -1 * maxTangentialVelocity * (1 + 0.5 / alpha ) * (coreRadius ./ r2) .* (1 - exp(-1 * alpha  .* r2.^2 / coreRadius^2));

% Set NaN velocities (which occur at zero radius) equal to zero
uTheta1(isnan(uTheta1)) = 0;
uTheta2(isnan(uTheta2)) = 0;

% Calculate vorticity for the first core
vort1 =      maxTangentialVelocity * coreRadius * (1 + 0.5 / alpha) * 2 * alpha / coreRadius^2 .* exp(-1 * alpha .* r1.^2 / coreRadius^2);
vort2 = -1 * maxTangentialVelocity * coreRadius * (1 + 0.5 / alpha) * 2 * alpha / coreRadius^2 .* exp(-1 * alpha .* r2.^2 / coreRadius^2);

% Add the vorticities
VORTICITY = vort1 + vort2;

% Convert polar velocites to cartesian for first core
u1 = -1 * uTheta1 .* sin(th1);
v1 =      uTheta1 .* cos(th1);

% Convert polar velocites to cartesian for second core
u2 = -1 * uTheta2 .* sin(th2);
v2 =      uTheta2 .* cos(th2);

% Add the velocities together to get the total field. Reshape into row
% vectors.
u = reshape(u1 + u2, nPoints, 1);
v = reshape(v1 + v2, nPoints, 1);

% Concatonate into a single vector.
U = [u; v];
   
end






