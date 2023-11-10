function [out] = modelMe(P, ST, Z)
% MODELME This function processes spatial and angular data related to an 
% animal's movement and computes various statistical measures to model the 
% animal's behavior. It involves thresholding, behavioral segmentation, 
% generation of ratemaps, optimization of model parameters, and calculation 
% of variance explained and modulation strength.
%
%   Inputs:
%     P:  A position vector containing timestamps and coordinates in the 
%         form [t, x1, y1, x2, y2] for two points, or [t, x, y] for a 
%         single point, where t is time in seconds.
%
%     ST: A spike times vector indicating the times at which neuronal 
%         spikes occurred, in seconds.
%
%     Z:  An angular variable in degrees, which could represent head 
%         direction or movement direction, etc.
%
%   Outputs:
%     out: A struct containing the following fields:
%          - out.model: Contains the modeled ratemap, fitting error, 
%            optimized parameters, and optimization metadata.
%          - out.data: Contains the original ratemap, raw and normalized 
%            ratemaps, occupancy counts and times, percentage of bins 
%            preserved after speed and location thresholding.
%          - out.MVLmeanHD: Mean Vector Length for head direction tuning.
%          - out.MVLmeanRH: Mean Vector Length for reference-heading model 
%            fit tuning.
%
%   Detailed Description:
%     The function takes the positional data and spike times to perform 
%     several computational steps:
%     1. Thresholding based on speed to filter out periods of low mobility.
%     2. Segmentation of behavioral data into spatial and angular bins.
%     3. Generation of ratemaps that represent the firing rates across 
%        spatial and angular bins.
%     4. Optimization of model parameters to fit the ratemap data.
%     5. Computation of Modulation Strength through Mean Vector Length 
%        (MVL) analysis of head direction and model tuning curves.
%
%   Example Usage:
%     P = [timestamp, x_position, y_position];
%     ST = [spike_time1, spike_time2, ...];
%     Z = [angular_data_in_degrees];
%     result = modelMe(P, ST, Z);
%
%
% References: 
%   - Jercog et al., 2019
%
% Note: 
%   - Code assumes that angular data Z is initially within 0-360 degrees 
%     range before being shifted to -180 to 180 degree range.
%   - Speed is calculated using an external function `get_speed`.
%   - Optimization initial conditions can be randomly chosen or set to the 
%     same values as in Jercog et al., 2019.

% Set Thresholds
thresh_spd = 4; % cm/s - speed threshold for movement

% Format inputs
t = P(:,1); % Extract time from the position vector
x = P(:,2); % Extract x-coordinates from the position vector
y = P(:,3); % Extract y-coordinates from the position vector
tpf = mode(diff(t)); % Calculate the most common time difference between frames to estimate time per frame

% Define start & end times
startTime = t(1); 
stopTime = t(end);

% Remove spike times that are outside of the time range of the position data
ST = ST(ST < stopTime & ST > startTime);

% Shift angular data Z to range [-180, 180] degrees for standardization
Z(Z > 180) = Z(Z > 180) - 360;

% Create a histogram of spikes for each time bin
t_edges = linspace(startTime, stopTime, numel(t)+1)';
SpkTrn = histcounts(ST, t_edges)';

% Calculate the speed using the external function `get_speed`
spd = get_speed(x, y, t);

% Identify indices where speed is below the threshold (out-of-bounds)
oob = spd < thresh_spd;

% Remove data where speed is below the threshold to exclude periods of low mobility
t(oob) = [];
x(oob) = [];
y(oob) = [];
SpkTrn(oob) = [];
Z(oob) = [];

% Segment behavioral data into spatial bins
nBins = 10; % Number of spatial bins
[~, xEdges, yEdges, binX, binY] = histcounts2(x, y, nBins);

% Calculate the center of each bin
xCenter = (diff(xEdges) / 2) + xEdges(1:end-1);
yCenter = (diff(yEdges) / 2) + yEdges(1:end-1);

% Segment angular data into bins
num_Z_bins = 10; % Number of angular bins
Z_bins = linspace(0, 360, num_Z_bins + 1);
Z_bin_ctrs = ((diff(Z_bins) / 2) + Z_bins(1:end-1)) - 180;


%% GENERATE RATEMAPS
% Preallocate data structures
r_xy = nan(nBins, num_Z_bins);
r_xyh = nan(nBins, nBins, num_Z_bins);
R_xyh = nan(nBins, nBins, num_Z_bins);
time_H = nan(nBins, nBins, num_Z_bins);
count_H = nan(nBins, nBins, num_Z_bins);
Occ = cell(nBins, nBins);
x_rc = zeros(nBins);
y_rc = zeros(nBins);

% Iterate over each spatial bin and calculate ratemaps
for rr = 1:nBins
    for cc = 1:nBins
        % Find frames in which animal occupied this spatial bin
        idx_here = find(binY == rr & binX == cc);
        
        % Calculate time spent in this bin
        time_in_bin = numel(idx_here) * tpf; 
        
        % Filter spikes and angular data for the current spatial bin
        spikes_here = SpkTrn(idx_here); 
        Z_here = Z(idx_here);
        
        % Compute occupancy and count for each angular bin
        [Z_count_here, ~] = histcounts(Z_here, Z_edges);
        [Z_occ_here, ~] = histcounts(Z_here, Z_bins);
        
        % Store occupancy time for each angular bin
        Occ{rr, cc} = Z_occ_here * tpf;
        
        % Compute ratemaps
        R_xyh(rr, cc, :) = Z_count_here / time_in_bin;
        time_H(rr, cc, :) = time_in_bin;
        count_H(rr, cc, :) = Z_count_here;
        
        % Check for sufficient occupancy before calculating rates
        if time_in_bin >= bin_threshold
            r_xy(rr, cc) = sum(spikes_here) / time_in_bin; % Rate for spatial bin
            r_xyh(rr, cc, :) = Z_count_here ./ Z_occ_here; % Rate for each head direction bin within the spatial bin
        end
        
        % Fill bin centers
        x_rc(rr, cc) = xCenter(cc);
        y_rc(rr, cc) = yCenter(rr);
    end
end

%% OPTIMIZATION PROCESS

% Define minimum firing rate threshold for a bin to be considered (Hz)
rCutOff = 0.5;

% Randomly choose initial conditions based on the number of bins
initialParams = choose_initial_conditions(nBins);

% Initial conditions for the 4 parameters to be optimized
pFitInitial = [initialParams.g; initialParams.thetaP; initialParams.xref; initialParams.yref]; 

% Set options for fminsearch optimization algorithm:
% 'Display' - controls display of output (turned off for cleaner output)
% 'TolX' - termination tolerance on the solution (smaller value for higher precision)
% 'TolFun' - termination tolerance on the function value (smaller value for higher precision)
% 'MaxFunEvals' - maximum number of function evaluations (increased for extensive search)
% 'MaxIter' - maximum number of iterations (increased for extensive search)
options = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8, ...
                   'MaxFunEvals', 8000, 'MaxIter', 8000);

% Perform the optimization using fminsearch by fitting the model to the data
% The anonymous function defines the cost to be minimized
[pFit, ~, exitflag, output] = fminsearch(@(pFit) FitLM(pFit, R_xyh, r_xy, rCutOff, nBins), ...
    pFitInitial, options);

% Retrieve model-predicted firing rates for the best-fit parameters
[R_xyh_model, fF] = get_Rxyh_model(pFit, R_xyh, r_xy, rCutOff, nBins);

% Match NaN values in the model predictions and empirical data to ensure
% comparability when calculating model performance or error
[RDatanan, rDatanan, rSDatanan] = matchnans(R_xyh_model, R_xyh, r_xyh);

%% MODULATION STRENGTH ANALYSIS

% Initialize matrices
% NaN indicates bins that don't exceed the firing rate cutoff or lack finite data
MVL = nan(nBins, nBins);
MVL_RH = nan(nBins, nBins);
mu = nan(nBins, nBins);
mu_RH = nan(nBins, nBins);
angCount = 0; angCount_rh = 0;

% Iterate over all spatial bins to calculate modulation strength
for rr = 1:nBins
    for cc = 1:nBins
        % Only proceed if average firing rate in the spatial bin exceeds cutoff
        if (rSDatanan(rr, cc) > rCutOff)
            
            % Retrieve the angular tuning curve for the current spatial bin
            rT = squeeze(RDatanan(rr, cc, :));
            rT_rh = squeeze(R_xyh_model(rr, cc, :));
            
            % Find indices of finite values in both head direction and RH tuning curves
            iF = find(isfinite(rT));
            iF_rh = find(isfinite(rT_rh));
            
            % Process head direction tuning curve if it contains finite values
            if ~isempty(iF)
                % Extract finite values from the tuning curve and corresponding angular bins
                tc_now = rT(iF);
                bins_now = deg2rad(Z_bin_ctrs(iF))';  % Convert from degrees to radians
                
                % Calculate the circular mean and mean resultant length for
                % HD (this is the "data" tuning curve)
                mu(rr, cc) = rad2deg(circ_mean(bins_now, tc_now));
                MVL(rr, cc) = circ_r(bins_now, tc_now);
                
                % Update counter for bins exceeding cutoff
                angCount = angCount + numel(iF);
            end
            
            % Process reference-heading tuning curve if it contains finite values
            if ~isempty(iF_rh)
                % Extract finite values from the RH tuning curve and corresponding angular bins
                tc_now_rh = rT_rh(iF_rh);
                bins_now_rh = deg2rad(Z_bin_ctrs(iF_rh))';  % Convert from degrees to radians
                
                % Calculate the circular mean and mean resultant length for
                % RH (this is the "model" tuning curve)
                mu_RH(rr, cc) = rad2deg(circ_mean(bins_now_rh, tc_now_rh));
                MVL_RH(rr, cc) = circ_r(bins_now_rh, tc_now_rh);
                
                % Update counter for RH bins exceeding cutoff
                angCount_rh = angCount_rh + numel(iF_rh);
            end
        end
    end
end

% Calculate the linear average of tuning strengths across spatial bins
% This represents the neuron's 'tuning strength' for the entire environment
MVL(MVL==Inf) = NaN; 
MVLmeanHD = nanmean(MVL, 'all');

MVL_RH(MVL_RH==Inf) = NaN;
MVLmeanRH = nanmean(MVL_RH, 'all');


%% PREPARE OUTPUTS
% MODEL CLASS
out.model.Rxyh = R_xyh_model;
out.model.error = fF;
out.model.fitParams.g = pFit(1);
out.model.fitParams.thetaP = pFit(2);
out.model.fitParams.xref = pFit(3);
out.model.fitParams.yref = pFit(4);
out.model.optim.exitflag = exitflag;
out.model.optim.output = output;
out.model.optim.options = options;

% DATA CLASS
out.data.Rxyh = R_xyh;
out.data.rxyh = r_xyh;
out.data.rxy = r_xy;
out.data.RxyhN = RDatanan;
out.data.rxyhN = rDatanan;
out.data.rxyN = rSDatanan;
out.data.occ.count = count_H;
out.data.occ.countmat = Occ;
out.data.occ.time = time_H;

% MEASURES
out.measures.MVL.HD = MVL;
out.measures.MVL.RH = MVL_RH;
out.measures.TS.HD = MVLmeanHD;
out.measures.TS.RH = MVLmeanRH;
out.measures.mu.HDCorr = MVL;
out.measures.mu.RHCorr = MVL_RH;

% INFO
out.info.bin.X = x_rc;
out.info.bin.Y = y_rc;

end