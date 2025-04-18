%============== INITIAL CONDITIONS =============%
function [Sf,Sp,Mf,Mp,Lp] = sub_init_fish_hist_noD(ID,DAYS,Sf,Sp,Mf,Mp,Lp)

%===== VARIABLES =====%
    %%%! Number of spatial cells
    NX = length(ID);

    %%%! Fraction of time in pelagic
    Sf.td = ones(NX,1);
    Sp.td = ones(NX,1);
    Mf.td = ones(NX,1);
    Mp.td = ones(NX,1);
    Lp.td = ones(NX,1);

    %%%! Rates, etc.
    nzero = {'met' 'enc_f' 'enc_p' 'enc_d' 'enc_zm' 'enc_zl' 'enc_be' ...
        'con_f' 'con_p' 'con_d' 'con_zm' 'con_zl' 'con_be' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'egg' 'rep' 'rec' 'clev' 'caught'};
    for i = 1 : length(nzero)
        Sf.(nzero{i}) = zeros(NX,1);
        Sp.(nzero{i}) = zeros(NX,1);
        Mf.(nzero{i}) = zeros(NX,1);
        Mp.(nzero{i}) = zeros(NX,1);
        Lp.(nzero{i}) = zeros(NX,1);
    end

    %%%! Proportion of population spawning (needed for phenology runs)
    Sf.S = ones(NX,DAYS);
    Sp.S = ones(NX,DAYS);
    Mf.S = ones(NX,DAYS);
    Mp.S = ones(NX,DAYS);
    Lp.S = ones(NX,DAYS);


end
