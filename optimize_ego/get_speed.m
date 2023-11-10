function speed = get_speed(x, y, t)
% GET_SPEED calculates the instantaneous speed from positional and time data.
%
% This function computes the speed of an object moving in a 2D space.
% The speed at each time step is calculated using the change in position over time.
% The first speed value is duplicated to match the length of the input vectors.
%
% INPUTS:
%   x : A vector of x-coordinates (in units consistent with y).
%   y : A vector of y-coordinates (in units consistent with x).
%   t : A vector of time stamps (in seconds), same length as x and y.
%
% OUTPUT:
%   speed : A vector of instantaneous speeds (in units per second), same length as x and y.

% Calculate the differences in position.
d_x = diff(x);
d_y = diff(y);

% Calculate the mode of the time differences to use as a consistent time step.
d_t = mode(diff(t));

% Compute the Euclidean distance between consecutive positions.
d_distance = sqrt((d_x.^2) + (d_y.^2));

% Calculate the raw speed values by dividing the distance by the time step.
speed_raw = abs(d_distance ./ d_t);

% Deal with boundaries of the speed vector
speed = zeros(length(speed_raw) + 1, 1);
speed(2:end) = speed_raw;
speed(1) = NaN;

end




