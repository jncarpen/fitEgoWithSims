function initial = choose_initial_conditions(total_bins)
    % CHOOSE_INITIAL_CONDITIONS selects random initial conditions for the parameters
    % that are input to the RH-model (see: get_Rxyh_model.m).
    % 
    % INPUT:
    % total_bins      Number of bins for the discretization of the space.
    % 
    % OUTPUT:
    % initial         A struct with the following fields containing initial conditions:
    %  - .g           A random gain factor between 0 and 1.
    %  - .thetaP      A randomly chosen orientation from a specified range.
    %  - .xref        A randomly chosen x-coordinate from the possible positions.
    %  - .yref        A randomly chosen y-coordinate from the possible positions.

    % Define the range of possible positions for x and y, and orientations.
    x_bins = 1:0.5:total_bins;
    y_bins = 1:0.5:total_bins;
    orientation = -180:1:180;
    
    % Randomly sample one element from each of the possible positions and orientation.
    % The rand function is used to generate a random number for the gain factor,
    % whereas randsample is used to randomly select from the bins and orientation.
    howMany = 1;
    initial.g = rand(1); 
    initial.thetaP = randsample(orientation, howMany)'; 
    initial.xref = randsample(x_bins, howMany)';       
    initial.yref = randsample(y_bins, howMany)';
end
