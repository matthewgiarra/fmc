function TUCKEYWINDOW = tuckey2(FILTERDIMENSIONS, SCALING)

% Create a hann (i.e. cosine) window
hannWindow = hann2(FILTERDIMENSIONS);

% Scale the window
TUCKEYWINDOW = SCALING * hannWindow;

% Cut off the values greater than 1
TUCKEYWINDOW(TUCKEYWINDOW > 1) = 1;

end