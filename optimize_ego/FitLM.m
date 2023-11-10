function [fF] = FitLM(pFit,rF,rP,rCutOff,Nbins)
% FitLM Computes the fit of a model using mean squared error.
%
% This function is a wrapper around get_Rxyh_model and it is used during the
% optimization process to evaluate the goodness of fit of the model to the data.
% It only outputs the mean squared error component of what get_Rxyh_model returns.
%
% INPUTS:
%   pFit     : A vector of initial parameters for the model fitting.
%   rF       : A 3D matrix of firing rates as a function of spatial bin and head direction.
%   rP       : A 2D matrix of average firing rates for each spatial bin.
%   rCutOff  : A scalar threshold for the minimum firing rate to consider a bin.
%   Nbins    : The number of bins used for discretizing the space.
%
% OUTPUT:
%   fF       : The mean squared error of the fit.

% Call the model fitting function and discard its first output using ~
[~, fF] = get_Rxyh_model(pFit, rF, rP, rCutOff, Nbins);
end
