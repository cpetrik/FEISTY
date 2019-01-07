%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,ENVR] = sub_futbio(ID,DY,NEMURO,ENVR,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,dfrate)

global DAYS GRD NX
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit h gam
global bent_eff rfrac Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel
global tstep K CGRD ni nj

%%% If biomass < individual fish mass per grid cell, set all rates to zero? %%%

%%% NEMURO information
ENVR = get_NEMURO(NEMURO,ID,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);
ENVR.dZm = sub_neg(ENVR.dZm);
ENVR.dZl = sub_neg(ENVR.dZl);

% Update benthic biomass with new detritus avail at that time step
[BENT.mass,BENT.pred] = sub_update_be(ENVR.Tb,BENT.mass,bent_eff,ENVR.det,[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);
BENT.mass = sub_check(BENT.mass);

% Pelagic-demersal coupling
%Lp: fraction of time large piscivores spends in pelagic
%Ld: fraction of time large demersals spends in pelagic
if (pdc == 0)
    Lp.td = ones(NX,1);
    Ld.td = zeros(NX,1);
elseif (pdc == 1)
    Lp.td = ones(NX,1);
    Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (pdc == 2)
    Lp.td = sub_tdif_pel(ENVR.H,Mf.bio,Mp.bio,Md.bio);
    Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
else
    Lp.td = ones(NX,1);
    Ld.td = ones(NX,1);
end

% Metabolism
Sf.met = sub_met(ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.met = sub_met(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Sd.met = sub_met(ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Mf.met = sub_met(ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.met = sub_met(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Md.met = sub_met(ENVR.Tp,ENVR.Tb,Md.td,M_m);
Lp.met = sub_met(ENVR.Tp,ENVR.Tb,Lp.td,M_l);
Ld.met = sub_met(ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Encounter rates
%           sub_enc(Tp     ,Tb     ,wgt,prey   ,tpel ,tprey,pref)
Sf.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sf.td,Sf.td,1);
Sp.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sp.td,Sp.td,1);
Sd.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sd.td,Sd.td,1);

Mf.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mf.td,Mf.td,MF_phi_MZ);
Mf.enc_zl = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mf.td,Mf.td,MF_phi_LZ);
Mf.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sf.bio,Mf.td,Mf.td,MF_phi_S);
Mf.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mf.td,Mf.td,MF_phi_S);
Mf.enc_d  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sd.bio,Mf.td,Mf.td,MF_phi_S);

Mp.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mp.td,Mp.td,MP_phi_MZ);
Mp.enc_zl = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mp.td,Mp.td,MP_phi_LZ);
Mp.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sf.bio,Mp.td,Mp.td,MP_phi_S);
Mp.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mp.td,Mp.td,MP_phi_S);
Mp.enc_d  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sd.bio,Mp.td,Mp.td,MP_phi_S);

Md.enc_be = sub_enc(ENVR.Tp,ENVR.Tb,M_m,BENT.mass,Md.td,1-Md.td,MD_phi_BE);

Lp.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mf.bio,Lp.td,Lp.td,LP_phi_MF);
Lp.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mp.bio,Lp.td,Lp.td,LP_phi_MP);
Lp.enc_d  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Md.bio,Lp.td,1-Lp.td,LP_phi_MD);

Ld.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mf.bio,Ld.td,Ld.td,LD_phi_MF);
Ld.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mp.bio,Ld.td,Ld.td,LD_phi_MP);
Ld.enc_d  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Md.bio,Ld.td,1-Ld.td,LD_phi_MD);
Ld.enc_be = sub_enc(ENVR.Tp,ENVR.Tb,M_l,BENT.mass,Ld.td,1-Ld.td,LD_phi_BE);

% Consumption rates
Sf.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sf.td,M_s,Sf.enc_zm);
Sp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sp.td,M_s,Sp.enc_zm);
Sd.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sd.td,M_s,Sd.enc_zm);

Mf.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_zl,Mf.enc_zm,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_f,Mf.enc_zm,Mf.enc_zl,Mf.enc_p,Mf.enc_d]);
Mf.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_p,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_d]);
Mf.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_d,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p]);

Mp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zl,Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_f,Mp.enc_zm,Mp.enc_zl,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_d]);
Mp.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_d,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p]);

Md.con_be = sub_cons(ENVR.Tp,ENVR.Tb,Md.td,M_m,Md.enc_be);

Lp.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_f,Lp.enc_p,Lp.enc_d]);
Lp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);
Lp.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_d,Lp.enc_p,Lp.enc_f]);

Ld.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_f,Ld.enc_p,Ld.enc_d,Ld.enc_be]);
Ld.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_p,Ld.enc_f,Ld.enc_d,Ld.enc_be]);
Ld.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_d,Ld.enc_p,Ld.enc_f,Ld.enc_be]);
Ld.con_be = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_be,Ld.enc_f,Ld.enc_p,Ld.enc_d]);

% Offline coupling
%MZ consumption cannot exceed amount lost to higher predation in NEMURO runs
[Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
    sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);
%LZ consumption cannot exceed amount lost to higher predation in NEMURO runs
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
Sf.clev = sub_clev(Sf.I,ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.clev = sub_clev(Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Sd.clev = sub_clev(Sd.I,ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Mf.clev = sub_clev(Mf.I,ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.clev = sub_clev(Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Md.clev = sub_clev(Md.I,ENVR.Tp,ENVR.Tb,Md.td,M_m);
Lp.clev = sub_clev(Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,M_l);
Ld.clev = sub_clev(Ld.I,ENVR.Tp,ENVR.Tb,Ld.td,M_l);

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
Sf.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Sd.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Mf.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Md.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Md.td,M_m);
Lp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Lp.td,M_l);
Ld.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Energy available for somatic growth nu
[Sf.nu, Sf.prod] = sub_nu(Sf.I,Sf.bio,Sf.met);
[Sp.nu, Sp.prod] = sub_nu(Sp.I,Sp.bio,Sp.met);
[Sd.nu, Sd.prod] = sub_nu(Sd.I,Sd.bio,Sd.met);
[Mf.nu, Mf.prod] = sub_nu(Mf.I,Mf.bio,Mf.met);
[Mp.nu, Mp.prod] = sub_nu(Mp.I,Mp.bio,Mp.met);
[Md.nu, Md.prod] = sub_nu(Md.I,Md.bio,Md.met);
[Lp.nu, Lp.prod] = sub_nu(Lp.I,Lp.bio,Lp.met);
[Ld.nu, Ld.prod] = sub_nu(Ld.I,Ld.bio,Ld.met);

% Maturation (note subscript on Kappa is larvae, juv, adult)
Sf.gamma = sub_gamma(K_l,Z_s,Sf.nu,Sf.die,Sf.bio,Sf.nmort,0,0);
Sp.gamma = sub_gamma(K_l,Z_s,Sp.nu,Sp.die,Sp.bio,Sp.nmort,0,0);
Sd.gamma = sub_gamma(K_l,Z_s,Sd.nu,Sd.die,Sd.bio,Sd.nmort,0,0);
Mf.gamma = sub_gamma(K_a,Z_m,Mf.nu,Mf.die,Mf.bio,Mf.nmort,dfrate,MFsel);
Mp.gamma = sub_gamma(K_j,Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,dfrate,MPsel);
Md.gamma = sub_gamma(K_j,Z_m,Md.nu,Md.die,Md.bio,Md.nmort,dfrate,MDsel);
Lp.gamma = sub_gamma(K_a,Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,dfrate,LPsel);
Ld.gamma = sub_gamma(K_a,Z_l,Ld.nu,Ld.die,Ld.bio,Ld.nmort,dfrate,LDsel);

% Egg production (by med and large size classes only)
[Mf.gamma,Mf.nu,Mf.rep,Mf.egg] = sub_rep(Mf.gamma,Mf.nu,K_a,Mf.S(:,DY),Mf.egg);
[Lp.gamma,Lp.nu,Lp.rep,Lp.egg] = sub_rep(Lp.gamma,Lp.nu,K_a,Lp.S(:,DY),Lp.egg);
[Ld.gamma,Ld.nu,Ld.rep,Ld.egg] = sub_rep(Ld.gamma,Ld.nu,K_a,Ld.S(:,DY),Ld.egg);

% Recruitment (from smaller size class)
Sf.rec = sub_rec_larv(Mf.rep,Mf.bio,rfrac);
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,rfrac);
Sd.rec = sub_rec_larv(Ld.rep,Ld.bio,rfrac);
Mf.rec = sub_rec(Sf.gamma,Sf.bio);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Md.rec = sub_rec(Sd.gamma,Sd.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);
Ld.rec = sub_rec(Md.gamma,Md.bio);

% Mass balance
Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.egg,Sf.nmort);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.egg,Sp.nmort);
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.egg,Sd.nmort);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.egg,Mf.nmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.egg,Mp.nmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.egg,Md.nmort);

Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.egg,Lp.nmort);
Ld.bio = sub_update_fi(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.egg,Ld.nmort);

% Fishing by rate
[Mf.bio, Mf.caught, Mf.fmort] = sub_fishing_rate(Mf.bio,dfrate,MFsel);
[Mp.bio, Mp.caught, Mp.fmort] = sub_fishing_rate(Mp.bio,dfrate,MPsel);
[Md.bio, Md.caught, Md.fmort] = sub_fishing_rate(Md.bio,dfrate,MDsel);
[Lp.bio, Lp.caught, Lp.fmort] = sub_fishing_rate(Lp.bio,dfrate,LPsel);
[Ld.bio, Ld.caught, Ld.fmort] = sub_fishing_rate(Ld.bio,dfrate,LDsel);
% [Mf.bio, Mf.caught, Mf.fmort] = sub_quad_fishing(Mf.bio,dfrate,MFsel,ENVR.Tp,ENVR.Tb,Mf.td);
% [Mp.bio, Mp.caught, Mp.fmort] = sub_quad_fishing(Mp.bio,dfrate,MPsel,ENVR.Tp,ENVR.Tb,Mp.td);
% [Md.bio, Md.caught, Md.fmort] = sub_quad_fishing(Md.bio,dfrate,MDsel,ENVR.Tp,ENVR.Tb,Md.td);
% [Lp.bio, Lp.caught, Lp.fmort] = sub_quad_fishing(Lp.bio,dfrate,LPsel,ENVR.Tp,ENVR.Tb,Lp.td);
% [Ld.bio, Ld.caught, Ld.fmort] = sub_quad_fishing(Ld.bio,dfrate,LDsel,ENVR.Tp,ENVR.Tb,Ld.td);

% Advection-Diffusion
% Sf.bio = sub_advec_vel(CGRD,Sf.bio,ENVR.U,ENVR.V,ni,nj,tstep);
% Sp.bio = sub_advec_vel(CGRD,Sp.bio,ENVR.U,ENVR.V,ni,nj,tstep);
% Sd.bio = sub_advec_vel(CGRD,Sd.bio,ENVR.U,ENVR.V,ni,nj,tstep);
% Mf.bio = sub_diff_sep(CGRD,Mf.bio,K,ni,nj,tstep);
% Mp.bio = sub_diff_sep(CGRD,Mp.bio,K,ni,nj,tstep);
% Md.bio = sub_diff_sep(CGRD,Md.bio,K,ni,nj,tstep);
% Lp.bio = sub_diff_sep(CGRD,Lp.bio,K,ni,nj,tstep);
% Ld.bio = sub_diff_sep(CGRD,Ld.bio,K,ni,nj,tstep);
% [Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,Md.bio,Lp.bio,Ld.bio] = sub_diff(CGRD,K,...
%     ni,nj,tstep,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,Md.bio,Lp.bio,Ld.bio);

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
