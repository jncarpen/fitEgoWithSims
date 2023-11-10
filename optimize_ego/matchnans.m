function [RDatanan, rDatanan, rSDatanan] = matchnans(R_xyh_model, R_xyh, r_xyh)
% MATCHNANS Propagates NaN values from a model matrix to corresponding data matrices.
%
% This function takes the NaN values from a 3D model matrix and applies the same NaN
% positions to two corresponding data matrices. It also calculates the average non-NaN 
% values for each [row, col] across the third dimension.
%
% INPUTS:
%   R_xyh_model : A 3D matrix representing the model, with potential NaN values.
%   R_xyh       : A 3D matrix representing observed data corresponding to R_xyh_model.
%   r_xyh       : Another 3D matrix representing observed data corresponding to R_xyh_model.
%
% OUTPUTS:
%   RDatanan    : The R_xyh matrix with NaN values matched to those in R_xyh_model.
%   rDatanan    : The r_xyh matrix with NaN values matched to those in R_xyh_model.
%   rSDatanan   : A 2D matrix with the average of non-NaN values from rDatanan for each [row, col].

% Get the size of the input matrices.
[r, c, v] = size(R_xyh_model);

% Flatten the 3D matrices to 1D for NaN-indexing.
R_xyh_model_flat = reshape(R_xyh_model, [], 1);
R_xyh_flat = reshape(R_xyh, [], 1);
r_xyh_flat = reshape(r_xyh, [], 1);

% Find indices of NaN in the model and apply NaN to those indices in the data matrices.
nanidx = isnan(R_xyh_model_flat);
R_xyh_flat(nanidx) = NaN;
r_xyh_flat(nanidx) = NaN;

% Reshape the flattened matrices back to their original 3D shape.
RDatanan = reshape(R_xyh_flat, r, c, v);
rDatanan = reshape(r_xyh_flat, r, c, v);

% Calculate the average non-NaN values.
rSDatanan = squeeze(nanmean(rDatanan, 3));
end