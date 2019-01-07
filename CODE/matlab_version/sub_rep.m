%%% BIOMASS MADE FROM REPRODUCTION
function [gamma, nu, rep, egg] = sub_rep(gamma,nu,K,S,egg)
%nu: energy for growth or spawning
%K: proportion allocated to growth
%S: fraction = fraction of pop spawning at that time
%S: 0 or 1 = indicates if spawning season (1=yes)
%egg: energy stored for later repro

global NX

% NOTE: Still never going to accumulate biomass as muscle tissue
% If it is spawning season, it gets spawned
% If it is not spawning season, it gets stored as repro energy "egg"

    if K<1.0
        rho = zeros(NX,1);
        id = (nu > 0.0);
        rho(id,1) = (1.0-K) .* nu(id,1);  %energy available for from eating
        
        id2 = (S>0.0);
        %rep = eggs from energy now + eggs from stored energy
        rep(id2,1) = S(id2,1) .* (rho(id2,1)+egg(id2,1));         %fraction of pop reproducing now
        egg(id2,1) = (1.0-S(id2,1)) .* (rho(id2,1)+egg(id2,1));   %rest gets stored for later
        
        id3 = (S<=0.0);
        rep(id3,1) = 0.0;
        egg(id3,1) = egg(id3,1) + rho(id3,1);
        
        %add what would be growth to next size up as repro
        rep = rep + gamma;
        gamma = zeros(NX,1);
        
        id4 = (nu > 0.0);
        %nu is now split into used for repro (nu) and stored (egg)
        nu(id4,1) = rep(id4,1);
    else
        rep = zeros(NX,1);
        egg = zeros(NX,1);
    end


end
