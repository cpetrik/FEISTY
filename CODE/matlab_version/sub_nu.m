%%% ENERGY AVAILABLE FOR GROWTH NU
function [nu, prod] = sub_nu(I,B,met)

    global Lambda
    
    % convert to biomass specific ingestion
    %nu = ((I/B)*Lambda) - met
    %nu = 0.5*I
    % Already in biomass specific ingestion
    nu = (I*Lambda) - met;
    prod = nu .* B;
end
