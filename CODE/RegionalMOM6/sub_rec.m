%%% Biomass recruiting to size-class (g m-2 d-1)
function rec = sub_rec(X,bio)
    % X = biomass specific maturation rate of smaller size class (gamma)
    % bio = biomass of smaller size class
    rec = X .* bio;
end
