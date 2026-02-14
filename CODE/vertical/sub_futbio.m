%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,ENVR] = sub_futbio(DY,ESM,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,param)

dfrate = param.dfrate;

%%% ESM information
ENVR = get_ESM(ESM,param,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);
ENVR.dZm  = sub_neg(ENVR.dZm);
ENVR.dZl  = sub_neg(ENVR.dZl);

if isfield(ENVR,'H')==0
    ENVR.H = param.depth;
end

if isfield(ENVR,'dz')==0
    ENVR.dz = param.dz;
end

%ENVR.det = ENVR.det *10;

% Update benthic biomass with new detritus avail at that time step
[BENT.mass,BENT.pred] = sub_update_be(BENT.mass,param,ENVR.det,[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);
BENT.mass = sub_check(BENT.mass);

% Pelagic-demersal coupling
%depth levels occupied
if (param.pdc == 1)
    Md.td = zeros(param.NX,1);
    Ld.td = zeros(param.NX,1);
    Md.td(param.NX) = 1.0;
    Ld.td(param.NX) = 1.0;
    if (ENVR.H <= param.PI_be_cutoff)
        Ld.td = ones(param.NX,1);
    end
end

% Metabolism
Sf.met = sub_met(ENVR.Tp,Sf.td,param.M_s,param);
Sp.met = sub_met(ENVR.Tp,Sp.td,param.M_s,param);
Sd.met = sub_met(ENVR.Tp,Sd.td,param.M_s,param);
Mf.met = sub_met(ENVR.Tp,Mf.td,param.M_m,param);
Mp.met = sub_met(ENVR.Tp,Mp.td,param.M_m,param);
Md.met = sub_met(ENVR.Tp,Md.td,param.M_m,param);
Lp.met = sub_met(ENVR.Tp,Lp.td,param.M_l,param);
Ld.met = sub_met(ENVR.Tp,Ld.td,param.M_l,param);

% Encounter rates
%           sub_enc(params,Tp     ,wgt     ,prey   ,tpel ,pref)
Sf.enc_zm = sub_enc(param,ENVR.Tp,param.M_s,ENVR.Zm,Sf.td,1);
Sp.enc_zm = sub_enc(param,ENVR.Tp,param.M_s,ENVR.Zm,Sp.td,1);
Sd.enc_zm = sub_enc(param,ENVR.Tp,param.M_s,ENVR.Zm,Sd.td,1);

Mf.enc_zm = sub_enc(param,ENVR.Tp,param.M_m,ENVR.Zm,Mf.td,param.MF_phi_MZ);
Mf.enc_zl = sub_enc(param,ENVR.Tp,param.M_m,ENVR.Zl,Mf.td,param.MF_phi_LZ);
Mf.enc_f  = sub_enc(param,ENVR.Tp,param.M_m,Sf.bio,Mf.td,param.MF_phi_S);
Mf.enc_p  = sub_enc(param,ENVR.Tp,param.M_m,Sp.bio,Mf.td,param.MF_phi_S);
Mf.enc_d  = sub_enc(param,ENVR.Tp,param.M_m,Sd.bio,Mf.td,param.MF_phi_S);

Mp.enc_zm = sub_enc(param,ENVR.Tp,param.M_m,ENVR.Zm,Mp.td,param.MP_phi_MZ);
Mp.enc_zl = sub_enc(param,ENVR.Tp,param.M_m,ENVR.Zl,Mp.td,param.MP_phi_LZ);
Mp.enc_f  = sub_enc(param,ENVR.Tp,param.M_m,Sf.bio,Mp.td,param.MP_phi_S);
Mp.enc_p  = sub_enc(param,ENVR.Tp,param.M_m,Sp.bio,Mp.td,param.MP_phi_S);
Mp.enc_d  = sub_enc(param,ENVR.Tp,param.M_m,Sd.bio,Mp.td,param.MP_phi_S);

Md.enc_be = sub_enc(param,ENVR.Tp,param.M_m,BENT.mass,Md.td,param.MD_phi_BE/100);

Lp.enc_f  = sub_enc(param,ENVR.Tp,param.M_l,Mf.bio,Lp.td,param.LP_phi_MF);
Lp.enc_p  = sub_enc(param,ENVR.Tp,param.M_l,Mp.bio,Lp.td,param.LP_phi_MP);
Lp.enc_d  = sub_enc(param,ENVR.Tp,param.M_l,Md.bio,Lp.td,param.LP_phi_MD);

Ld.enc_f  = sub_enc(param,ENVR.Tp,param.M_l,Mf.bio,Ld.td,param.LD_phi_MF);
Ld.enc_p  = sub_enc(param,ENVR.Tp,param.M_l,Mp.bio,Ld.td,param.LD_phi_MP);
Ld.enc_d  = sub_enc(param,ENVR.Tp,param.M_l,Md.bio,Ld.td,param.LD_phi_MD/100);
Ld.enc_be = sub_enc(param,ENVR.Tp,param.M_l,BENT.mass,Ld.td,param.LD_phi_BE/100);

% Consumption rates
Sf.con_zm = sub_cons(param,ENVR.Tp,Sf.td,param.M_s,Sf.enc_zm);
Sp.con_zm = sub_cons(param,ENVR.Tp,Sp.td,param.M_s,Sp.enc_zm);
Sd.con_zm = sub_cons(param,ENVR.Tp,Sd.td,param.M_s,Sd.enc_zm);

Mf.con_zm = sub_cons(param,ENVR.Tp,Mf.td,param.M_m,[Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_zl = sub_cons(param,ENVR.Tp,Mf.td,param.M_m,[Mf.enc_zl,Mf.enc_zm,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_f  = sub_cons(param,ENVR.Tp,Mf.td,param.M_m,[Mf.enc_f,Mf.enc_zm,Mf.enc_zl,Mf.enc_p,Mf.enc_d]);
Mf.con_p  = sub_cons(param,ENVR.Tp,Mf.td,param.M_m,[Mf.enc_p,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_d]);
Mf.con_d  = sub_cons(param,ENVR.Tp,Mf.td,param.M_m,[Mf.enc_d,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p]);

Mp.con_zm = sub_cons(param,ENVR.Tp,Mp.td,param.M_m,[Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_zl = sub_cons(param,ENVR.Tp,Mp.td,param.M_m,[Mp.enc_zl,Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_f  = sub_cons(param,ENVR.Tp,Mp.td,param.M_m,[Mp.enc_f,Mp.enc_zm,Mp.enc_zl,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(param,ENVR.Tp,Mp.td,param.M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_d]);
Mp.con_d  = sub_cons(param,ENVR.Tp,Mp.td,param.M_m,[Mp.enc_d,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p]);

Md.con_be = sub_cons(param,ENVR.Tp,Md.td,param.M_m,Md.enc_be);

Lp.con_f  = sub_cons(param,ENVR.Tp,Lp.td,param.M_l,[Lp.enc_f,Lp.enc_p,Lp.enc_d]);
Lp.con_p  = sub_cons(param,ENVR.Tp,Lp.td,param.M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);
Lp.con_d  = sub_cons(param,ENVR.Tp,Lp.td,param.M_l,[Lp.enc_d,Lp.enc_p,Lp.enc_f]);

Ld.con_f  = sub_cons(param,ENVR.Tp,Ld.td,param.M_l,[Ld.enc_f,Ld.enc_p,Ld.enc_d,Ld.enc_be]);
Ld.con_p  = sub_cons(param,ENVR.Tp,Ld.td,param.M_l,[Ld.enc_p,Ld.enc_f,Ld.enc_d,Ld.enc_be]);
Ld.con_d  = sub_cons(param,ENVR.Tp,Ld.td,param.M_l,[Ld.enc_d,Ld.enc_p,Ld.enc_f,Ld.enc_be]);
Ld.con_be = sub_cons(param,ENVR.Tp,Ld.td,param.M_l,[Ld.enc_be,Ld.enc_f,Ld.enc_p,Ld.enc_d]);

% Offline coupling
if param.hploss == 1
    %MZ consumption cannot exceed amount lost to higher predation in COBALT runs
    [Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
        sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);
    %LZ consumption cannot exceed amount lost to higher predation in COBALT runs
    [Mf.con_zl,Mp.con_zl,ENVR.fZl] = ...
        sub_offline_zl(Mf.con_zl,Mp.con_zl,Mf.bio,Mp.bio,ENVR.dZl);
end
%Track fraction of benthic material consumed
[ENVR.fB] = sub_offline_bent([Md.con_be,Ld.con_be],[Md.bio,Ld.bio],BENT.mass);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sf.I = Sf.con_zm;
Sp.I = Sp.con_zm;
Sd.I = Sd.con_zm;
Mf.I = Mf.con_zm + Mf.con_zl + Mf.con_f + Mf.con_p + Mf.con_d;
Mp.I = Mp.con_zm + Mp.con_zl + Mp.con_f + Mp.con_p + Mp.con_d;
Md.I = Md.con_be;
Lp.I = Lp.con_f + Lp.con_p + Lp.con_d;
Ld.I = Ld.con_f + Ld.con_p + Ld.con_d + Ld.con_be;

% Consumption related to Cmax
Sf.clev = sub_clev(param,Sf.I,ENVR.Tp,Sf.td,param.M_s);
Sp.clev = sub_clev(param,Sp.I,ENVR.Tp,Sp.td,param.M_s);
Sd.clev = sub_clev(param,Sd.I,ENVR.Tp,Sd.td,param.M_s);
Mf.clev = sub_clev(param,Mf.I,ENVR.Tp,Mf.td,param.M_m);
Mp.clev = sub_clev(param,Mp.I,ENVR.Tp,Mp.td,param.M_m);
Md.clev = sub_clev(param,Md.I,ENVR.Tp,Md.td,param.M_m);
Lp.clev = sub_clev(param,Lp.I,ENVR.Tp,Lp.td,param.M_l);
Ld.clev = sub_clev(param,Ld.I,ENVR.Tp,Ld.td,param.M_l);

% Death rates (g m-2 d-1)
Sf.die = Mp.con_f.*Mp.bio + Mf.con_f.*Mf.bio;
Sp.die = Mp.con_p.*Mp.bio + Mf.con_p.*Mf.bio;
Sd.die = Mp.con_d.*Mp.bio + Mf.con_d.*Mf.bio;
Mf.die = Lp.con_f.*Lp.bio + Ld.con_f.*Ld.bio;
Mp.die = Lp.con_p.*Lp.bio + Ld.con_p.*Ld.bio;
Md.die = Lp.con_d.*Lp.bio + Ld.con_d.*Ld.bio;

% predation rates (m-2 d-1)
Sf.pred = Sf.die ./ Sf.bio;
Sp.pred = Sp.die ./ Sp.bio;
Sd.pred = Sd.die ./ Sd.bio;
Mf.pred = Mf.die ./ Mf.bio;
Mp.pred = Mp.die ./ Mp.bio;
Md.pred = Md.die ./ Md.bio;

% Natural mortality rates
Sf.nmort = sub_nmort(param,ENVR.Tp,Sf.td,param.M_s);
Sp.nmort = sub_nmort(param,ENVR.Tp,Sp.td,param.M_s);
Sd.nmort = sub_nmort(param,ENVR.Tp,Sd.td,param.M_s);
Mf.nmort = sub_nmort(param,ENVR.Tp,Mf.td,param.M_m);
Mp.nmort = sub_nmort(param,ENVR.Tp,Mp.td,param.M_m);
Md.nmort = sub_nmort(param,ENVR.Tp,Md.td,param.M_m);
Lp.nmort = sub_nmort(param,ENVR.Tp,Lp.td,param.M_l);
Ld.nmort = sub_nmort(param,ENVR.Tp,Ld.td,param.M_l);

% Energy available for somatic growth nu
[Sf.nu, Sf.prod] = sub_nu(param,Sf.I,Sf.bio,Sf.met);
[Sp.nu, Sp.prod] = sub_nu(param,Sp.I,Sp.bio,Sp.met);
[Sd.nu, Sd.prod] = sub_nu(param,Sd.I,Sd.bio,Sd.met);
[Mf.nu, Mf.prod] = sub_nu(param,Mf.I,Mf.bio,Mf.met);
[Mp.nu, Mp.prod] = sub_nu(param,Mp.I,Mp.bio,Mp.met);
[Md.nu, Md.prod] = sub_nu(param,Md.I,Md.bio,Md.met);
[Lp.nu, Lp.prod] = sub_nu(param,Lp.I,Lp.bio,Lp.met);
[Ld.nu, Ld.prod] = sub_nu(param,Ld.I,Ld.bio,Ld.met);


% Maturation (note subscript on Kappa is larvae, juv, adult)
Sf.gamma = sub_gamma(param.K_l,param.Z_s,Sf.nu,Sf.die,Sf.bio,Sf.nmort,0,0);
Sp.gamma = sub_gamma(param.K_l,param.Z_s,Sp.nu,Sp.die,Sp.bio,Sp.nmort,0,0);
Sd.gamma = sub_gamma(param.K_l,param.Z_s,Sd.nu,Sd.die,Sd.bio,Sd.nmort,0,0);
Mf.gamma = sub_gamma(param.K_a,param.Z_m,Mf.nu,Mf.die,Mf.bio,Mf.nmort,dfrate,param.MFsel);
Mp.gamma = sub_gamma(param.K_j,param.Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,dfrate,param.MPsel);
Md.gamma = sub_gamma(param.K_j,param.Z_m,Md.nu,Md.die,Md.bio,Md.nmort,dfrate,param.MDsel);
Lp.gamma = sub_gamma(param.K_a,param.Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,dfrate,param.LPsel);
Ld.gamma = sub_gamma(param.K_a,param.Z_l,Ld.nu,Ld.die,Ld.bio,Ld.nmort,dfrate,param.LDsel);


% Egg production (by med and large size classes only)
[Mf.gamma,Mf.nu,Mf.rep] = sub_rep(param.NX,Mf.gamma,Mf.nu,param.K_a);
[Lp.gamma,Lp.nu,Lp.rep] = sub_rep(param.NX,Lp.gamma,Lp.nu,param.K_a);
[Ld.gamma,Ld.nu,Ld.rep] = sub_rep(param.NX,Ld.gamma,Ld.nu,param.K_a);


% Recruitment (from smaller size class)
Sf.rec = sub_rec_larv(Mf.rep,Mf.bio,param.rfrac);
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,param.rfrac);

Mf.rec = sub_rec(Sf.gamma,Sf.bio);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);

Sd.rec = sub_rec_larvSD(Ld.rep,Ld.bio,param.rfrac,ENVR);
Md.rec = sub_recMD(Sd.gamma,Sd.bio,Md.td,ENVR);
Ld.rec = sub_recLD(Md.gamma,Md.bio,ENVR);

% Fishing by rate
[Mf.caught, Mf.fmort] = sub_fishing_rate(Mf.bio,dfrate,param.MFsel);
[Mp.caught, Mp.fmort] = sub_fishing_rate(Mp.bio,dfrate,param.MPsel);
[Md.caught, Md.fmort] = sub_fishing_rate(Md.bio,dfrate,param.MDsel);
[Lp.caught, Lp.fmort] = sub_fishing_rate(Lp.bio,dfrate,param.LPsel);
[Ld.caught, Ld.fmort] = sub_fishing_rate(Ld.bio,dfrate,param.LDsel);

% Mass balance
%                      (bio_in,  rec,   nu,   rep,   gamma,   die,   nmort,fmort)
Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.nmort,0);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.nmort,0);
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.nmort,0);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.nmort,Mf.fmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.nmort,Mp.fmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.nmort,Md.fmort);

Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.nmort,Lp.fmort);
Ld.bio = sub_update_fi(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.nmort,Ld.fmort);


% Forward Euler checks for demographics and movement
Sf.bio=sub_check(Sf.bio);
Sp.bio=sub_check(Sp.bio);
Sd.bio=sub_check(Sd.bio);
Mf.bio=sub_check(Mf.bio);
Mp.bio=sub_check(Mp.bio);
Md.bio=sub_check(Md.bio);
Lp.bio=sub_check(Lp.bio);
Ld.bio=sub_check(Ld.bio);

end
