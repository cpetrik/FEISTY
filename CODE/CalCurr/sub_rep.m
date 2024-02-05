%%% BIOMASS MADE FROM REPRODUCTION
function [gamma, nu, rep] = sub_rep(NX,gamma,nu,K)
%nu: energy for growth or spawning
%K: proportion allocated to growth

% NOTE: Still never going to accumulate biomass as muscle tissue

    if K<1.0
        rho = zeros(NX,1);
        id = (nu > 0.0);
        rho(id,1) = (1.0-K) .* nu(id,1);  %energy available for from eating
        
        %add what would be growth to next size up as repro
        rep = rho + gamma;
        gamma = zeros(NX,1);
        
    else
        rep = zeros(NX,1);
    end


end
