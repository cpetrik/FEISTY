%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sf,Sm,Sp,Sl,Sd,Mf,Mm,Mp,Ml,Md,Lp,Ll,Ld,BENT,ENVR] = sub_futbio_midw(ID,DY,COBALT,GRD,Sf,Sm,Sp,Sl,Sd,Mf,Mm,Mp,Ml,Md,Lp,Ll,Ld,BENT,param)

%%% If biomass < individual fish mass per grid cell, set all rates to zero? %%%

%%% COBALT information
ENVR = get_COBALT_midw(COBALT,GRD,ID,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);
ENVR.dZm = sub_neg(ENVR.dZm);
ENVR.dZl = sub_neg(ENVR.dZl);

% Update benthic biomass with new detritus avail at that time step
[BENT.mass,BENT.pred] = sub_update_be(BENT.mass,param,ENVR.det,[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);
BENT.mass = sub_check(BENT.mass);

% Pelagic-demersal coupling
%td = fraction of time spent in epipelagic, mesopelagic, benthic
if (param.pdc == 0)
    Lp.td = ones(param.NX,1);
    Ld.td = zeros(param.NX,1);
elseif (param.pdc == 1)
    Lp.td = ones(param.NX,1);
    Ld.td = sub_tdif_dem(ENVR.H,param,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (param.pdc == 2)
    Lp.td = sub_tdif_pel(ENVR.H,param,Mf.bio,Mp.bio,Md.bio);
    Ld.td = sub_tdif_dem(ENVR.H,param,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (param.pdc == 3)
    LP.td = sub_tdif_pmid(ENVR.H,param,Mf.bio,Mm.bio,Mp.bio,Ml.bio,Md.bio,BENT.mass);
    LL.td = sub_tdif_lmid(ENVR.H,param,Mf.bio,Mm.bio,Mp.bio,Ml.bio,Md.bio,BENT.mass);
    LD.td = sub_tdif_dmid(ENVR.H,param,Mf.bio,Mm.bio,Mp.bio,Ml.bio,Md.bio,BENT.mass);
else
    Lp.td = ones(param.NX,1);
    Ld.td = ones(param.NX,1);
end

% Metabolism
Sf.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Sf.td,param.M_s,param);
Sp.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Sp.td,param.M_s,param);
Sm.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Sm.td,param.M_s,param);
Sl.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Sl.td,param.M_s,param);
Sd.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Sd.td,param.M_s,param);
Mf.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Mf.td,param.M_m,param);
Mp.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Mp.td,param.M_m,param);
Mm.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Mm.td,param.M_m,param);
Ml.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Ml.td,param.M_m,param);
Md.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Md.td,param.M_m,param);
Lp.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Lp.td,param.M_l,param);
Ll.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Ll.td,param.M_l,param);
Ld.met = sub_met_midw(ENVR.Tp,ENVR.Tm,ENVR.Tb,Ld.td,param.M_l,param);

% Encounter rates
% Epipel overlap with zooplankton (assume 50-50 day-night)
%          Day            +  Night
zepi = 0.5*(1-param.zmigr) + 0.5*1.0;
% Mesopel overlap with zooplankton
zmeso = 0.5*(param.zmigr) + 0.5*1.0;
%           sub_enc(Depth,params,Tp     ,Tm     ,Tb     ,wgt      ,prey   ,tdif ,tprey,pref)
Sf.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_s,ENVR.Zm,Sf.td,zepi,1.0);
Sp.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_s,ENVR.Zm,Sp.td,zepi,1.0);
Sm.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_s,ENVR.Zm,Sm.td,zmeso,1.0);
Sl.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_s,ENVR.Zm,Sl.td,zmeso,1.0);
Sd.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_s,ENVR.Zm,Sd.td,zepi,1.0);

Mf.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zm,Mf.td,zepi,param.MF_phi_MZ);
Mf.enc_zl = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zl,Mf.td,zepi,param.MF_phi_LZ);
Mf.enc_f  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sf.bio,Mf.td,1.0,param.MF_phi_S);
Mf.enc_p  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sp.bio,Mf.td,1.0,param.MF_phi_S);
Mf.enc_d  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sd.bio,Mf.td,1.0,param.MF_phi_S);
Mf.enc_m  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sm.bio,Mf.td,1.0,param.MF_phi_S);
Mf.enc_l  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sl.bio,Mf.td,1.0,param.MF_phi_S);

Mp.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zm,Mp.td,zepi,param.MP_phi_MZ);
Mp.enc_zl = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zl,Mp.td,zepi,param.MP_phi_LZ);
Mp.enc_f  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sf.bio,Mp.td,1.0,param.MP_phi_S);
Mp.enc_p  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sp.bio,Mp.td,1.0,param.MP_phi_S);
Mp.enc_d  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sd.bio,Mp.td,1.0,param.MP_phi_S);
Mp.enc_m  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sm.bio,Mp.td,1.0,param.MP_phi_S);
Mp.enc_l  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sl.bio,Mp.td,1.0,param.MP_phi_S);

Mm.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zm,Mm.td,zmeso,param.MM_phi_MZ);
Mm.enc_zl = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zl,Mm.td,zmeso,param.MM_phi_LZ);
Mm.enc_f  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sf.bio,Mm.td,1.0,param.MM_phi_S);
Mm.enc_p  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sp.bio,Mm.td,1.0,param.MM_phi_S);
Mm.enc_d  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sd.bio,Mm.td,1.0,param.MM_phi_S);
Mm.enc_m  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sm.bio,Mm.td,1.0,param.MM_phi_S);
Mm.enc_l  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sl.bio,Mm.td,1.0,param.MM_phi_S);

Ml.enc_zm = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zm,Mm.td,zmeso,param.ML_phi_MZ);
Ml.enc_zl = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,ENVR.Zl,Mm.td,zmeso,param.ML_phi_LZ);
Ml.enc_f  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sf.bio,Ml.td,1.0,param.ML_phi_S);
Ml.enc_p  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sp.bio,Ml.td,1.0,param.ML_phi_S);
Ml.enc_d  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sd.bio,Ml.td,1.0,param.ML_phi_S);
Ml.enc_m  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sm.bio,Ml.td,1.0,param.ML_phi_S);
Ml.enc_l  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,Sl.bio,Ml.td,1.0,param.ML_phi_S);

Md.enc_be = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_m,BENT.mass,Md.td,1-Md.td,param.MD_phi_BE);

Lp.enc_f  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,Mf.bio,Lp.td,Lp.td,param.LP_phi_MF);
Lp.enc_p  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,Mp.bio,Lp.td,Lp.td,param.LP_phi_MP);
Lp.enc_d  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,Md.bio,Lp.td,1-Lp.td,param.LP_phi_MD);

Ld.enc_f  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,Mf.bio,Ld.td,Ld.td,param.LD_phi_MF);
Ld.enc_p  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,Mp.bio,Ld.td,Ld.td,param.LD_phi_MP);
Ld.enc_d  = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,Md.bio,Ld.td,1-Ld.td,param.LD_phi_MD);
Ld.enc_be = sub_enc(ENVR.H,param,ENVR.Tp,ENVR.Tm,ENVR.Tb,param.M_l,BENT.mass,Ld.td,1-Ld.td,param.LD_phi_BE);

% Consumption rates
Sf.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s,Sf.enc_zm);
Sp.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s,Sp.enc_zm);
Sd.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s,Sd.enc_zm);

Mf.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_zl = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_zl,Mf.enc_zm,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_f,Mf.enc_zm,Mf.enc_zl,Mf.enc_p,Mf.enc_d]);
Mf.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_p,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_d]);
Mf.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_d,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p]);

Mp.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_zl = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_zl,Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_f,Mp.enc_zm,Mp.enc_zl,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_d]);
Mp.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_d,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p]);

Md.con_be = sub_cons(param,ENVR.Tp,ENVR.Tb,Md.td,param.M_m,Md.enc_be);

Lp.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,[Lp.enc_f,Lp.enc_p,Lp.enc_d]);
Lp.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);
Lp.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,[Lp.enc_d,Lp.enc_p,Lp.enc_f]);

Ld.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_f,Ld.enc_p,Ld.enc_d,Ld.enc_be]);
Ld.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_p,Ld.enc_f,Ld.enc_d,Ld.enc_be]);
Ld.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_d,Ld.enc_p,Ld.enc_f,Ld.enc_be]);
Ld.con_be = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_be,Ld.enc_f,Ld.enc_p,Ld.enc_d]);

% Offline coupling
%MZ consumption cannot exceed amount lost to higher predation in COBALT runs
[Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
    sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);
%LZ consumption cannot exceed amount lost to higher predation in COBALT runs
[Mf.con_zl,Mp.con_zl,ENVR.fZl] = ...
    sub_offline_zl(Mf.con_zl,Mp.con_zl,Mf.bio,Mp.bio,ENVR.dZl);
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
Sf.clev = sub_clev(param,Sf.I,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s);
Sp.clev = sub_clev(param,Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s);
Sd.clev = sub_clev(param,Sd.I,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s);
Mf.clev = sub_clev(param,Mf.I,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m);
Mp.clev = sub_clev(param,Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m);
Md.clev = sub_clev(param,Md.I,ENVR.Tp,ENVR.Tb,Md.td,param.M_m);
Lp.clev = sub_clev(param,Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l);
Ld.clev = sub_clev(param,Ld.I,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l);

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
Sf.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s);
Sp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s);
Sd.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s);
Mf.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m);
Mp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m);
Md.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Md.td,param.M_m);
Lp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l);
Ld.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l);

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
Mf.gamma = sub_gamma(param.K_a,param.Z_m,Mf.nu,Mf.die,Mf.bio,Mf.nmort,param.dfrate,param.MFsel);
Mp.gamma = sub_gamma(param.K_j,param.Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,param.dfrate,param.MPsel);
Md.gamma = sub_gamma(param.K_j,param.Z_m,Md.nu,Md.die,Md.bio,Md.nmort,param.dfrate,param.MDsel);
Lp.gamma = sub_gamma(param.K_a,param.Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,param.dfrate,param.LPsel);
Ld.gamma = sub_gamma(param.K_a,param.Z_l,Ld.nu,Ld.die,Ld.bio,Ld.nmort,param.dfrate,param.LDsel);

% Egg production (by med and large size classes only)
[Mf.gamma,Mf.nu,Mf.rep] = sub_rep(param.NX,Mf.gamma,Mf.nu,param.K_a);
[Lp.gamma,Lp.nu,Lp.rep] = sub_rep(param.NX,Lp.gamma,Lp.nu,param.K_a);
[Ld.gamma,Ld.nu,Ld.rep] = sub_rep(param.NX,Ld.gamma,Ld.nu,param.K_a);

% Recruitment (from smaller size class)
Sf.rec = sub_rec_larv(Mf.rep,Mf.bio,param.rfrac);
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,param.rfrac);
Sd.rec = sub_rec_larv(Ld.rep,Ld.bio,param.rfrac);
Mf.rec = sub_rec(Sf.gamma,Sf.bio);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Md.rec = sub_rec(Sd.gamma,Sd.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);
Ld.rec = sub_rec(Md.gamma,Md.bio);

% Mass balance
Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.nmort);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.nmort);
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.nmort);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.nmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.nmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.nmort);

Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.nmort);
Ld.bio = sub_update_fi(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.nmort);

% Fishing by rate
[Mf.bio, Mf.caught, Mf.fmort] = sub_fishing_rate(Mf.bio,param.dfrate,param.MFsel);
[Mp.bio, Mp.caught, Mp.fmort] = sub_fishing_rate(Mp.bio,param.dfrate,param.MPsel);
[Md.bio, Md.caught, Md.fmort] = sub_fishing_rate(Md.bio,param.dfrate,param.MDsel);
[Lp.bio, Lp.caught, Lp.fmort] = sub_fishing_rate(Lp.bio,param.dfrate,param.LPsel);
[Ld.bio, Ld.caught, Ld.fmort] = sub_fishing_rate(Ld.bio,param.dfrate,param.LDsel);

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
