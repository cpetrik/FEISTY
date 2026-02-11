%============== INITIAL CONDITIONS =============%
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT] = sub_init_fish_multiyr(NX,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT)

%===== VARIABLES =====%

%%%! Fish biomass
    % Sf.bio = repmat(Sf.bio,NX,1) ./ NX;
    % Sp.bio = repmat(Sp.bio,NX,1) ./ NX;
    % Sd.bio = repmat(Sd.bio,NX,1) ./ NX;
    % Mf.bio = repmat(Mf.bio,NX,1) ./ NX;
    % Mp.bio = repmat(Mp.bio,NX,1) ./ NX;
    % Lp.bio = repmat(Lp.bio,NX,1) ./ NX;
    % 
    % Md.bio = repmat(Md.bio,NX,1);
    % Ld.bio = repmat(Ld.bio,NX,1);
    % 
    % %Maybe not keep track of whole water column, just seafloor?
    % Md.bio(1:(NX-1)) = zeros;
    % %Needs an if statement for seafloor depth at location
    % if param.depth <= 200
    %     Ld.bio = Ld.bio ./ NX;
    % else
    %     Ld.bio(1:(NX-1)) = zeros;
    % end

    %%%! Fraction of time in pelagic
    Sf.td = ones(NX,1);
    Sp.td = ones(NX,1);
    Sd.td = ones(NX,1);
    Mf.td = ones(NX,1);
    Mp.td = ones(NX,1);
    Md.td = zeros(NX,1);
    Lp.td = ones(NX,1);
    Ld.td = zeros(NX,1);
    Md.td(NX) = 1.0;
    Ld.td(NX) = 1.0;

    %%%! Rates, etc.
    nzero = {'met' 'enc_f' 'enc_p' 'enc_d' 'enc_zm' 'enc_zl' 'enc_be' ...
        'con_f' 'con_p' 'con_d' 'con_zm' 'con_zl' 'con_be' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'egg' 'rep' 'rec' 'clev' 'caught' 'fmort' ...
        'mort'};
    for i = 1 : length(nzero)
        Sf.(nzero{i}) = zeros(NX,1);
        Sp.(nzero{i}) = zeros(NX,1);
        Sd.(nzero{i}) = zeros(NX,1);
        Mf.(nzero{i}) = zeros(NX,1);
        Mp.(nzero{i}) = zeros(NX,1);
        Md.(nzero{i}) = zeros(NX,1);
        Lp.(nzero{i}) = zeros(NX,1);
        Ld.(nzero{i}) = zeros(NX,1);
    end
    
    %%%! Proportion of population spawning (needed for phenology runs)
    % Sf.S = ones(NX,DAYS);
    % Sp.S = ones(NX,DAYS);
    % Sd.S = ones(NX,DAYS);
    % Mf.S = ones(NX,DAYS);
    % Mp.S = ones(NX,DAYS);
    % Md.S = ones(NX,DAYS);
    % Lp.S = ones(NX,DAYS);
    % Ld.S = ones(NX,DAYS);

    %%%! Detritus 
    %maybe not keep track of whole water column, just seafloor?
    BENT.mass = zeros(NX,1);
    BENT.pred = zeros(NX,1);

    BENT.mass(NX) = BENT.bio;

end
