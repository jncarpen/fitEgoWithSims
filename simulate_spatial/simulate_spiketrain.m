function varargout = simulate_spiketrain(param, mode)
% SIMULATE_SPIKETRAIN Simulates a spiketrain for a neuron.
% This function generates spike times and spike train vectors based on the chosen
% simulation mode which could be spatial ('place'), egocentric head direction ('ego'), 
% or a combination of both ('placeego'). The output varies depending on the mode.
%
% Syntax:
%   [spiketimes, spiketrain] = simulate_spiketrain(param, 'place')
%   [spiketimes, spiketrain, pref_theta] = simulate_spiketrain(param, 'placeego')
%   [spiketimes, spiketrain] = simulate_spiketrain(param, 'ego')
%
% Inputs:
%   param - A structure containing simulation parameters.
%   mode  - A string specifying the simulation mode.
%
% Outputs (variable):
%   spiketimes  - A vector of simulated spike times.
%   spiketrain  - A binary vector indicating whether a spike occurred at each timestamp.
%   pref_theta  - Preferred direction offset (radians) (only for 'placeego' and 'ego' mode).
%
% Examples:
%   param = struct('map', rateMap, 'P', positionMatrix, 'noise', noiseLevel);
%   [times, train] = simulate_spiketrain(param, 'place');
%
%   param = struct('map', rateMap, 'P', positionMatrix, 'theta', 90, ...
%                  'Z', Directions, 'kappa', 5, 'noise', noiseLevel, ...
%                  'rp', referencePoint, 'A', peakRate);
%   [times, train, theta] = simulate_spiketrain(param, 'placeego');
%
% Notes:
%   - The 'place' mode generates spikes based on a spatial rate map only.
%   - The 'placeego' mode uses both the spatial rate map and direction.
%   - The 'ego' mode generates spikes based on direction only.
%   - Any circular variable can be input as param.Z (e.g., head or movement
%   direction)
%
% See also: simulate_place, simulate_placeego, simulate_ego, and
% simulate_ratemap (for calculation of the 'map' structure)

% Validate the mode and parameters
validateMode(mode);

% Call the appropriate external function based on the mode
switch lower(mode)
    case 'place'
        requiredFields = {'map', 'P', 'noise'};
        validateParameters(requiredFields);
        [spiketimes, spiketrain] = simulate_place(param.map, param.P, param.noise);
        varargout = {spiketimes, spiketrain};

    case 'placeego'
        requiredFields = {'map', 'P', 'noise', 'Z', 'A', 'theta', 'rp', 'kappa'};
        validateParameters(requiredFields); 
        [spiketimes, spiketrain, pref_theta] = simulate_placeego(param.map, param);
        varargout = {spiketimes, spiketrain, pref_theta};

    case 'ego'
        requiredFields = {'P', 'noise', 'Z', 'A', 'theta', 'rp', 'kappa'};
        validateParameters(requiredFields);
        [spiketimes, spiketrain, pref_theta] = simulate_ego(param);
        varargout = {spiketimes, spiketrain, pref_theta};

    otherwise
        error('Unsupported mode. Choose ''place'', ''placeego'', or ''ego''.');
end

% Helper function to validate the mode is a string
function validateMode(inputMode)
    if ~ischar(inputMode)
        error('Mode must be a string.');
    end
end

% Helper function to validate required parameters are present for the mode
function validateParameters(requiredFields)
    missingFields = setdiff(requiredFields, fieldnames(param));
    if ~isempty(missingFields)
        error(['Missing required parameters for mode ' mode ': ' strjoin(missingFields, ', ')]);
    end
end
end
