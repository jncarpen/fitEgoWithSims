function [spiketimes, spiketrain] = simulate_place(map, P, noise)
% SIMULATE_PLACE Simulates spike times based on a spatial rate map and head LEDs positions.
%
% This function generates spiketimes for a place cell model based on a specified 
% spatial rate map, a matrix P containing timestamps and head LED positions, and 
% a specified noise level. Spikes are generated using a Poisson process, with the rate
% determined by the rate map value corresponding to the animal's position at each timestamp.
%
% Inputs:
%   map        - A structure containing the rate map and bin occupancy data, or
%                a matrix representing the rate map (Hz) if not provided in a structure.
%   P          - An Nx5 matrix where the first column is timestamps (s) and the next 
%                four columns are the (x1, y1) and (x2, y2) positions of the head LEDs.
%   noise      - A scalar representing constant noise to be added to the rate.
%
% Outputs:
%   spiketimes - A vector of simulated spike times.
%   spiketrain - A binary vector indicating spike occurrence.
%
% Notes:
%   When specifying a place field center in the root structure as [a,b],
%   the center of the resulting place cell will be flipped, at [b,a].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the timestamp vector from the first column of P
t = P(:, 1);

% Calculate the sampling frequency based on the mode of the time differences
fs = mode(diff(t));

% Validate the rate map format and extract it
if isstruct(map); RM = map.z; else; RM = map; end

% Convert rate map to spikes per frame by multiplying by the sampling interval
lambdaMatrix = RM .* fs;

% Initialize the spiketrain and spiketimes vectors
spiketrain = zeros(length(t), 1);
spiketimes = [];

% Generate spikes for each timestamp
for frame = 1:length(t)
    
    % Determine the current bin based on the animal's position
    x_t = map.whichBin.x(frame);
    y_t = map.whichBin.y(frame);

    % Continue only if the current position is valid (i.e., not NaN)
    if ~isnan(x_t) && ~isnan(y_t)
        % Retrieve the rate map value for the current position and add noise
        lambda = lambdaMatrix(x_t, y_t) + noise;

        % Draw a random number of spikes from a Poisson distribution based on the rate
        n = poissrnd(lambda);

        % If spikes occurred, update the spiketrain and add the times to spiketimes
        if n > 0
            spiketrain(frame) = n;
            spiketimes = [spiketimes; repmat(t(frame), n, 1)];
        end
    end
end
end
