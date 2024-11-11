function generate_thresholds(sortedValues, num_iterations)
    % Iterate over the specified number of iterations
    for x = 1:num_iterations
        % Generate variable name dynamically with a suffix
        var_name = ['threshold_beta_' sprintf('%02d', x)];
        
        % Assign the threshold value to the dynamically generated variable
        threshold{s} = assignin('caller', var_name, sortedValues_beta(floor((x/10)*length(sortedValues_beta))));
        
    end
end