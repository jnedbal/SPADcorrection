function pointer = createPointer
% Function to create a camera focus pointer

%% Create a larger matrix of the pointer
pointer = NaN(32, 32);

%Work in polar coordinates.
% Create angular coordinates of one quadrant of the pointer
angle = 0 : pi / 192 : 2 * pi;
% Create radial coordinates of one quadrant of the pointer
radius = 13 : 0.1 : 15;
radius = reshape(radius, 1, 1, numel(radius));
radius = repmat(radius, [size(angle), 1]);
angle = repmat(angle, 1, 1, size(radius, 3));

% Create cartesian coordinates of one quadrant of the pointer
[x, y] = pol2cart(angle(:), radius(:));
% Round the cartesian coordinates and prevent rounding errors to zero,
% which would throw an error
x = round(x) + 16;
y = round(y) + 16;
x = max(x, 1);
y = max(y, 1);
% Convert all coordinates into a 2D map
for i = 1 : numel(x)
    pointer(x(i), y(i)) = 2;
end

% Work in polar coordinates.
% Create angular coordinates of one quadrant of the pointer
angle = ((-pi / 24) : pi / 192 : (pi / 24));
angle = angle + repmat((0 : pi / 4 : 15 * pi / 8)', 1, numel(angle));
% Create radial coordinates of one quadrant of the pointer
radius = 13 : 0.1 : 15;
radius = reshape(radius, 1, 1, numel(radius));
radius = repmat(radius, [size(angle), 1]);
angle = repmat(angle, 1, 1, size(radius, 3));

% Create cartesian coordinates of one quadrant of the pointer
[x, y] = pol2cart(angle(:), radius(:));
% Round the cartesian coordinates and prevent rounding errors to zero,
% which would throw an error
x = round(x) + 16;
y = round(y) + 16;
x = max(x, 1);
y = max(y, 1);
% Convert all coordinates into a 2D map
for i = 1 : numel(x)
    pointer(x(i), y(i)) = 1;
end

% combine the quadrants to make a circle
% pointer = [fliplr(pointer), pointer];
% pointer = [flipud(pointer); pointer];
%close all


% hard boundary
%figure(1)
%imagesc(pointer)
% weight the points: point itself; average of nearest neighbors;
% averaged of diagonal neighbors.  These must add up to 1.
% x = 0 : (size(pointer, 2) - 1);
% wp = .4;  wn = .4;  wd = .2;
% ind = 2:length(x)-1;
% pointer(ind,ind) = wp*pointer(ind,ind) ...
%   + (wn/4)*(pointer(ind-1,ind  ) + pointer(ind+1,ind  ) + pointer(ind  ,ind-1) + pointer(ind  ,ind+1) ) ...
%   + (wd/4)*(pointer(ind-1,ind-1) + pointer(ind-1,ind+1) + pointer(ind+1,ind-1) + pointer(ind+1,ind+1) );
% extended boundary
%figure(2)

% Get rid of the edges
%pointer = pointer(2 : end - 1, 2 : end - 1);

% % Get rid of the negative values
% pointer = pointer + 1;
% % Normalize to one
% pointer = pointer / 2;
% % Convert value 0.5 to NaN
%imagesc(pointer)