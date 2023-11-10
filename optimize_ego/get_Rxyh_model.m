function [RXYHM, fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins)
% get_Rxyh_model computes a model of firing rates across spatial and head direction bins
% and quantifies the goodness of this fit through mean squared error.
%
% INPUTS:
%   pFit:        Parameters for the model fit as a vector [g; thetaP; xref; yref], where:
%                g - Gain factor of the model.
%                thetaP - Preferred firing direction in degrees.
%                xref, yref - Reference point coordinates for angle calculations.
%
%   rF:          Three-dimensional matrix of firing rates R_data(x,y,H) with dimensions 
%                (Nbins x Nbins x Nbins), representing activity across spatial and head
%                direction bins.
%
%   rP:          Two-dimensional matrix of firing rates r(x,y) with dimensions 
%                (Nbins x Nbins), representing spatial bins.
%
%   rCutOff:     Minimum firing rate threshold for a bin to be included in the model.
%
%   Nbins:       Number of bins for each dimension (x, y, and H) used for discretization.
%
% OUTPUTS:
%   RXYHM:       Three-dimensional matrix (Nbins x Nbins x Nbins) containing the model's 
%                predicted firing rates for each spatial and head direction bin.
%
%   fF:          Mean squared error (MSE) between the model's predicted firing rates and 
%                the actual data, averaged over the bins that passed the firing rate 
%                threshold and had finite values.
% EXAMPLE USAGE:
%   [RXYHM, fF] = get_Rxyh_model(pFit, rF, rP, rCutOff, Nbins);
% References:
%   [1] Jercog et al.

% Extract parameters from input vector
g = pFit(1);      % Gain factor for model fit
thetaP = pFit(2); % Preferred direction in degrees (deg)
xref = pFit(3);   % X-coordinate of the reference point
yref = pFit(4);   % Y-coordinate of the reference point

% Initialize RXYHM matrix
RXYHM = nan(Nbins, Nbins, Nbins);

% Initialize sum of squared errors (SSE) and angular bin counts to zero
fF = 0;
angCount = 0;

% Iterate over each spatial bin to calculate model firing rates and SSE
for ii = 1:Nbins
    for jj = 1:Nbins
        
        % Check if the average firing rate in the current spatial bin exceeds the cutoff threshold
        if rP(ii, jj) > rCutOff
            % Extract the tuning curve for the current spatial bin (ii, jj) across all head direction bins
            rT = squeeze(rF(ii, jj, :));
            
            % Find indices of head direction bins that have finite values (ignoring NaN or Inf)
            iF = find(isfinite(rT));
            
            % Proceed only if there are finite values to work with
            if ~isempty(iF)
                % Calculate the allocentric angle for each head direction bin relative to the reference point
                % Mapping the result onto the range [-180, 180] degrees
                a = 180 * atan2(yref - ii, xref - jj) / pi - ...
                    (-180 - 360 / (2 * Nbins) + iF * 360 / Nbins);
                
                % Compute the cosine of the angular difference between egocentric bearing and preferred direction
                cFac = cos(pi*(a-thetaP)/180);
                
                % Calculate the expected value (mean) of the cosine factors
                cBar = mean(cFac, 'omitnan');

                % Scale and shift the cosine function by gain factor g and subtract the expected value
                z = 1 + g * (cFac - cBar);
                
                % Apply a Heaviside step function to retain only positive values
                z(z < 0) = 0;

                % Accumulate the sum of squared errors between the model and actual firing rates
                fF = fF + sum((z - rT(iF)).^2, 'omitnan');
                
                % Increment the count of angular bins that exceeded the firing rate threshold
                angCount = angCount + numel(iF);

                % Store the modeled firing rates in the RXYHM matrix for the current spatial and head direction bins
                RXYHM(ii, jj, iF) = z;
            end
        end
    end
end

% Compute the mean SSE by dividing by the number of angular bins that passed the threshold
fF = fF / angCount;

end