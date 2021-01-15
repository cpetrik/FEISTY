%============== INITIAL CONDITIONS =============%
function [fish,BENT] = sub_init_fishvec(ID,param)

%===== VARIABLES =====%
    %%%! Number of spatial cells
    NX = length(ID);

    %%%! Fish biomass
    X = 1.0e-5; % very small amount
    fish.bio = zeros(NX,length(param.wc));
    fish.bio(:,param.ixFish) = ones(NX,length(param.ixFish)) * X;

    %%%! Fraction of time in pelagic
    fish.td = ones(NX,length(param.wc));
    fish.td(:,[3 10 11]) = zeros;
    
    %%%! Rates, etc.
    nzero = {'met' 'enc' 'con' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'rep' 'rec' 'clev' 'caught' 'fmort'};
    for i = 1 : length(nzero)
        fish.(nzero{i}) = zeros(NX,length(param.wc));
    end
    
    %%%! Detritus
    BENT.mass = ones(NX,1) * X;
    BENT.pred = zeros(NX,1);

end
