%============== INITIAL CONDITIONS =============%
function [Sd,Md,Ld,BENT] = sub_init_benthos_be2(ID,DAYS)

%===== VARIABLES =====%
    %%%! Number of spatial cells
    NX = length(ID);

    %%%! Fish biomass
    X = 1.0e-5; % very small amount
    Sd.bio = ones(NX,1) * X;
    Md.bio = ones(NX,1) * X;
    Ld.bio = ones(NX,1) * X;

    %%%! Fraction of time in pelagic
    Sd.td = ones(NX,1);
    Md.td = zeros(NX,1);
    Ld.td = zeros(NX,1);

    %%%! Rates, etc.
    nzero = {'met' 'enc_f' 'enc_p' 'enc_d' 'enc_zm' 'enc_zl' 'enc_be' ...
        'con_f' 'con_p' 'con_d' 'con_zm' 'con_zl' 'con_be' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'egg' 'rep' 'rec' 'clev' 'caught' 'fmort'};
    for i = 1 : length(nzero)
        Sd.(nzero{i}) = zeros(NX,1);
        Md.(nzero{i}) = zeros(NX,1);
        Ld.(nzero{i}) = zeros(NX,1);
    end
    
    %%%! Proportion of population spawning (needed for phenology runs)
    Sd.S = ones(NX,DAYS);
    Md.S = ones(NX,DAYS);
    Ld.S = ones(NX,DAYS);
    
    %%%! Detritus
    BENT.sm = ones(NX,1) * X;
    BENT.md = ones(NX,1) * X;
    BENT.predS = zeros(NX,1);
    BENT.predM = zeros(NX,1);

end
