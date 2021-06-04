%============== INITIAL CONDITIONS =============%
function [Sf,Sm,Sp,Sl,Sd,Mf,Mm,Mp,Ml,Md,Lp,Ll,Ld,BENT] = sub_init_fish_midw(ID)

%===== VARIABLES =====%
    %%%! Number of spatial cells
    NX = length(ID);

    %%%! Fish biomass
    X = 1.0e-5; % very small amount
    Sf.bio = ones(NX,1) * X;
    Sm.bio = ones(NX,1) * X;
    Sp.bio = ones(NX,1) * X;
    Sl.bio = ones(NX,1) * X;
    Sd.bio = ones(NX,1) * X;
    Mf.bio = ones(NX,1) * X;
    Mm.bio = ones(NX,1) * X;
    Mp.bio = ones(NX,1) * X;
    Ml.bio = ones(NX,1) * X;
    Md.bio = ones(NX,1) * X;
    Lp.bio = ones(NX,1) * X;
    Ll.bio = ones(NX,1) * X;
    Ld.bio = ones(NX,1) * X;

    %%%! Fraction of time in pelagic
    tep = [1 0 0];
    tm  = [0.5 0.5 0];
    tl  = [0 1 0];
    tdm = [0 0 1];
    Sf.td = repmat(tep,NX,1);
    Sm.td = repmat(tep,NX,1);
    Sp.td = repmat(tep,NX,1);
    Sl.td = repmat(tep,NX,1);
    Sd.td = repmat(tep,NX,1);
    Mf.td = repmat(tep,NX,1);
    Mm.td = repmat(tm,NX,1);
    Mp.td = repmat(tep,NX,1);
    Ml.td = repmat(tl,NX,1);
    Md.td = repmat(tdm,NX,1);
    Lp.td = repmat(tep,NX,1);
    Ll.td = repmat(tl,NX,1);
    Ld.td = repmat(tdm,NX,1);



    %%%! Rates, etc.
    nzero = {'met' 'enc_f' 'enc_p' 'enc_d' 'enc_zm' 'enc_zl' 'enc_be' ...
        'con_f' 'con_p' 'con_d' 'con_zm' 'con_zl' 'con_be' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'rep' 'rec' 'clev' 'caught' 'fmort'};
    for i = 1 : length(nzero)
        Sf.(nzero{i}) = zeros(NX,1);
        Sm.(nzero{i}) = zeros(NX,1);
        Sp.(nzero{i}) = zeros(NX,1);
        Sl.(nzero{i}) = zeros(NX,1);
        Sd.(nzero{i}) = zeros(NX,1);
        Mf.(nzero{i}) = zeros(NX,1);
        Mm.(nzero{i}) = zeros(NX,1);
        Mp.(nzero{i}) = zeros(NX,1);
        Ml.(nzero{i}) = zeros(NX,1);
        Md.(nzero{i}) = zeros(NX,1);
        Lp.(nzero{i}) = zeros(NX,1);
        Ll.(nzero{i}) = zeros(NX,1);
        Ld.(nzero{i}) = zeros(NX,1);
    end

    %%%! Detritus
    BENT.mass = ones(NX,1) * X;
    BENT.pred = zeros(NX,1);

end
