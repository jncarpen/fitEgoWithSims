function [spiketimes, spiketrain, pref_theta] = simulate_ego(param)
% SIMULATE_EGO Generates spike times based on an egocentric bearing model.
%
% This function simulates spiketimes using a von Mises distribution to model
% the egocentric bearing of an animal in a 2D space. The spiking activity is
% modulated by the animal's orientation and proximity to a reference point in 
% the environment, simulating a head-direction or egocentric cell.
%
% Inputs:
%   param - A structure with the following fields:
%       theta   - Angular orientation of the animal relative to a reference point (degrees).
%       Z       - Angular variable representing the head direction at each timestamp (degrees).
%       P       - Nx3 matrix where the first column is timestamps (s), and the next two columns
%                 are the (x, y) positions of the animal.
%       kappa   - Concentration parameter of the von Mises distribution.
%       rp      - Reference point coordinates in the form [x y].
%       A       - Peak firing rate of the unit (Hz).
%       noise   - Constant noise value to be added to the spike rate (Hz).
%       nspikes - The desired number of spikes to simulate.
%
% Outputs:
%   spiketimes - A vector of simulated spike times.
%   spiketrain - A binary vector indicating whether a spike occurred at each timestamp.
%   pref_theta  - Preferred direction offset (radians).
%
% Notes:
%   The function uses a Poisson process to generate spikes, where the rate is
%   determined by the product of the spatial and angular probabilities, scaled
%   by the peak firing rate and the mode of the time difference between position samples.

% Unpack position and time from the param structure
x = param.P(:,2);
y = param.P(:,3);
t = param.P(:,1); 
fs = mode(diff(t));  % Time per frame (s)

% Convert the preferred direction from degrees to radians
% offset by 90 degrees to align 'up' as 0 radians.
pref_theta = deg2rad(mod(param.theta - 90, 360)) - pi;

% Calculate the egocentric bearing relative to the reference point
ego = deg2rad(mod(atan2d(param.rp(2) - y, param.rp(1) - x) - param.Z + 180, 360) - 180);

% Bin the egocentric bearing angles
[~, edges, bin] = histcounts(ego, linspace(-pi, pi, 101)); % Circular histogram
ctrs = (diff(edges) / 2) + edges(1:end-1);
bin(bin == 0) = NaN; 

% Construct the von Mises distribution
[vm_pdf, ~] = circ_vmpdf(ctrs, pref_theta, param.kappa);

% Scale the von Mises distribution by the peak firing rate and time per frame
vm_pdf = (param.A .* vm_pdf) .* fs;

% Initialize the spiketrain vector
spiketrain = zeros(length(t), 1);
spiketimes = [];

% Generate spikes for each timestamp
for frame = 1:length(t)
    
    % Skip if the current frame has NaN position
    if isnan(bin(frame))
        continue;
    end

    % Compute the lambda (spike rate) for the current frame
    lambda = vm_pdf(bin(frame)) + param.noise;

    % Draw a random number of spikes from a Poisson distribution based on the rate
    n = poissrnd(lambda);

    % If spikes occurred, update the spiketrain and add the times to spiketimes
    if n > 0
        spiketrain(frame) = n;
        spiketimes = [spiketimes; repmat(t(frame), n, 1)];
    end
end
end


