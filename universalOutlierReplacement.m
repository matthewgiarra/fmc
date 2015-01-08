function VALUE_REPLACED = universalOutlierReplacement(VALUE, ISOUTLIER, STENCILRADIUS)
% This function performs universal outlier replacement on a
% rectilinear scalar field. 

% Calculate the height and width of the vector field.
[fieldHeight, fieldWidth] = size(VALUE);

% Stencil radius for this pass
r = STENCILRADIUS;

% Outside stencil element numbers
stencilElements = [(1 : (2 * r + 1) * r + r), (2 * r + 1) * r + r + 2 : (2 * r + 1)^2];

% Replicate the candidate replacement values
VALUE_REPLACED = VALUE;

% Loop over the whole grid
for m = 1 + r : fieldHeight - r
   for n = 1 + r : fieldWidth - r
       
       % Extract the scalar value its surroundings values at
       % each grid point.
       valueStencil = VALUE(m - r : m + r, n - r : n + r);

       % Stencil as columns
       valueCol = valueStencil(:);

       % Stencil elements excluding center point
       valueColOutside = valueCol(stencilElements);

       % If the value at grid point (m, n) was flagged as an outlier,
       % then replace its value with the median of the surrounding values
       if ISOUTLIER(m, n)
          VALUE_REPLACED(m, n) = median(valueColOutside); 
       end       
   end
end

end


