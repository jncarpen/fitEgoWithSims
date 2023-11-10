function [map] = simulate_ratemap(root)
% SIMULATE_RATEMAP Generates a Gaussian rate map for a place field.
%
% This function creates a 2D Gaussian rate map based on specified parameters
% for a simulated place field within a square arena. It computes the firing rate
% distribution and associates each position with a bin number.
%
% Inputs:
%   root    - A structure with fields describing the place field and arena:
%             A:       Peak firing rate (Hz).
%             ctr:     Center of the place field (x,y).
%             sigma:   Spread of the place field (x,y).
%             size:    Size of the square arena (one side, in cm).
%             bins:    Number of bins along each dimension.
%             P:       Position vector (timestamps, x, y), optional.
%
% Outputs:
%   map     - A structure with fields:
%             z:        Rate map (Hz) representing firing rates at each bin.
%             pdx:      Normalized probability density function of the rate map.
%             whichBin: Structure with x and y fields containing the bin number for each timestamp.
%
% Example:
%   root = struct('A', 10, 'ctr', [50, 50], 'sigma', [10, 10], 'size', 100, 'bins', 20);
%   map = simulate_ratemap(root);

% Extract place field parameters from the input structure
A = root.A;  
x0 = root.ctr(1); 
y0 = root.ctr(2);  
sigmaX = root.sigma(1);  
sigmaY = root.sigma(2); 

% Define the bin edges and centers for the arena
bins = root.bins;  
sz = root.size;  
xEdges = linspace(0, sz, bins + 1);
yEdges = linspace(0, sz, bins + 1);
[xCenter, yCenter] = meshgrid((diff(xEdges) / 2) + xEdges(1:end-1), ...
                              (diff(yEdges) / 2) + yEdges(1:end-1));

% Determine the bin number for each position sample
xpos = root.P(:, 2); ypos = root.P(:, 3);
xbin = discretize(xpos, xEdges);
ybin = discretize(ypos, yEdges);

% Generate the Gaussian rate map (Hz)
fxy = exp(-(((xCenter - x0).^2) / (2 * sigmaX^2) + ...
                ((yCenter - y0).^2) / (2 * sigmaY^2)));

% Scale the rate map by the peak firing rate
fxy_scaled = A * fxy;

% Package the output in a structure
map.z = fxy_scaled;        % Rate map (Hz)
map.pdx = fxy;             % Probability density function of the rate map
map.whichBin.x = xbin;     % Bin number for x positions
map.whichBin.y = ybin;     % Bin number for y positions
end

