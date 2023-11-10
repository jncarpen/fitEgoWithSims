function [spiketimes, spiketrain, pref_theta] = simulate_placeego(map, param)
% SIMULATE_PLACEEGO Simulates spiketimes based on place and egocentric head direction.
%
% This function generates spiketimes using a combination of a spatial rate map and
% an egocentric head direction model. Spikes are generated using a Poisson process
% where the rate is determined by the product of spatial and head direction probabilities,
% scaled by the peak firing rate and constant noise.
%
% Inputs:
%   map   - A structure containing the rate map 'z' and probability density 'pdx',
%           or a matrix representing the rate map (Hz) if not provided in a structure.
%   param - A structure with the following fields:
%           A: Peak firing rate of the unit (Hz).
%           P: Nx3 matrix of position vectors [time, x, y].
%           theta: Head direction (degrees, 0 is up).
%           Z: Head direction for each timestamp (degrees, 0 should be 'up').
%           rp: Reference point in the form [x, y].
%           kappa: Concentration parameter for von Mises distribution.
%           noise: Constant noise value (Hz).
%
% Outputs:
%   spiketimes  - A vector of simulated spiketimes.
%   spiketrain  - A binary vector indicating whether a spike occurred at each timestamp.
%   pref_theta  - Preferred direction offset (radians).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack position and head direction data from the param structure
x = param.P(:,2);
y = param.P(:,3);
t = param.P(:,1);

% Calculate time per frame from the mode of the timestamp differences
fs = mode(diff(t)); % Sampling frequency (Hz)

% Convert the preferred direction from degrees to radians
% offset by 90 degrees to align 'up' as 0 radians.
pref_theta = deg2rad(mod(param.theta + 90, 360)) - pi;

% Calculate egocentric bearing from the allocentric angle and head direction
ego = deg2rad(mod(atan2d(param.rp(2) - y, param.rp(1) - x) - param.Z + 180, 360) - 180);

% Bin the egocentric bearing angles
[~, edges, bin] = histcounts(ego, linspace(-pi, pi, 101));
ctrs = (diff(edges) / 2) + edges(1:end-1);
bin(bin == 0) = NaN; 


% Construct the von Mises distribution (unscaled)   
[vm_pdf, ~] = circ_vmpdf(ctrs, pref_theta, param.kappa);

% Initialize the spiketrain vector
spiketrain = zeros(length(t), 1);
spiketimes = [];

% Generate spikes for each timestep
% based on the joint probability of the subject being in a given
% spatial and angular bin
for frame = 1:length(t)

    % Find current bin indices for x, y, and head direction
    x_t = map.whichBin.x(frame);
    y_t = map.whichBin.y(frame);
    z_t = bin(frame);

    % Check for valid data and compute the lambda (spike rate) for the current frame
    if ~isnan(x_t) || ~isnan(y_t) || ~isnan(z_t)
        p_angular = vm_pdf(z_t);
        p_spatial = map.pdx(x_t, y_t);
        lambda = ((param.A .* (p_angular * p_spatial)) .* fs) + param.noise;

        % Draw a random number of spikes from a Poisson distribution based on the rate
        n = poissrnd(lambda);

        if n > 0
            spiketrain(frame) = n;
            spiketimes = [spiketimes; repmat(t(frame), n, 1)];
        end 
    end
end
end
