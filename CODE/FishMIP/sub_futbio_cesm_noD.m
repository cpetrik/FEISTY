%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sf,Sp,Mf,Mp,Lp,ENVR] = sub_futbio_cesm_noD(ID,DY,CESM,ENVR,Sf,Sp,Mf,Mp,Lp,dfrate)

global DAYS GRD NX
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit h gam
global rfrac Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S 
global LP_phi_MF LP_phi_MP LP_phi_MD 
global MFsel MPsel LPsel LDsel

%%% If biomass < individual fish mass per grid cell, set all rates to zero? %%%

%%% CESM information
ENVR = get_CESM(CESM,ID,DY);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);

% Pelagic-demersal coupling
%Lp: fraction of time large piscivores spends in pelagic
%Ld: fraction of time large demersals spends in pelagic
if (pdc == 0)
    Lp.td = ones(NX,1);
%     Ld.td = zeros(NX,1);
elseif (pdc == 1)
    Lp.td = ones(NX,1);
%     Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (pdc == 2)
    Lp.td = sub_tdif_pel(ENVR.H,Mf.bio,Mp.bio,Md.bio);
%     Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
else
    Lp.td = ones(NX,1);
%     Ld.td = ones(NX,1);
end

% Metabolism
Sf.met = sub_met(ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.met = sub_met(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Mf.met = sub_met(ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.met = sub_met(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Lp.met = sub_met(ENVR.Tp,ENVR.Tb,Lp.td,M_l);

% Encounter rates
%           sub_enc(Tp     ,Tb     ,wgt,prey   ,tpel ,tprey,pref)
Sf.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sf.td,Sf.td,1);
Sp.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sp.td,Sp.td,1);

Mf.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mf.td,Mf.td,MF_phi_MZ);
Mf.enc_zl = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mf.td,Mf.td,MF_phi_LZ);
Mf.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sf.bio,Mf.td,Mf.td,MF_phi_S);
Mf.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mf.td,Mf.td,MF_phi_S);

Mp.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mp.td,Mp.td,MP_phi_MZ);
Mp.enc_zl = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mp.td,Mp.td,MP_phi_LZ);
Mp.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sf.bio,Mp.td,Mp.td,MP_phi_S);
Mp.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mp.td,Mp.td,MP_phi_S);

Lp.enc_f  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mf.bio,Lp.td,Lp.td,LP_phi_MF);
Lp.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mp.bio,Lp.td,Lp.td,LP_phi_MP);

% Consumption rates
Sf.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sf.td,M_s,Sf.enc_zm);
Sp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sp.td,M_s,Sp.enc_zm);

Mf.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_zl,Mf.enc_zm,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_f,Mf.enc_zm,Mf.enc_zl,Mf.enc_p,Mf.enc_d]);
Mf.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_p,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_d]);

Mp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zl,Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_f,Mp.enc_zm,Mp.enc_zl,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_d]);

Lp.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_f,Lp.enc_p,Lp.enc_d]);
Lp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);

% Offline coupling
% %MZ consumption cannot exceed amount lost to higher predation in CESM runs
% [Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
%     sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);
% %LZ consumption cannot exceed amount lost to higher predation in CESM runs
% [Mf.con_zl,Mp.con_zl,ENVR.fZl] = ...
%     sub_offline_zl(Mf.con_zl,Mp.con_zl,Mf.bio,Mp.bio,ENVR.dZl);
% %Track fraction of benthic material consumed
% [ENVR.fB] = sub_offline_bent([Md.con_be,Ld.con_be],[Md.bio,Ld.bio],BENT.mass);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sf.I = Sf.con_zm;
Sp.I = Sp.con_zm;
Mf.I = Mf.con_zm + Mf.con_zl + Mf.con_f + Mf.con_p;
Mp.I = Mp.con_zm + Mp.con_zl + Mp.con_f + Mp.con_p;
Lp.I = Lp.con_f + Lp.con_p;

% Consumption related to Cmax
Sf.clev = sub_clev(Sf.I,ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.clev = sub_clev(Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Mf.clev = sub_clev(Mf.I,ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.clev = sub_clev(Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Lp.clev = sub_clev(Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,M_l);

% Death rates (g m-2 d-1)
Sf.die = Mp.con_f.*Mp.bio + Mf.con_f.*Mf.bio;
Sp.die = Mp.con_p.*Mp.bio + Mf.con_p.*Mf.bio;
Mf.die = Lp.con_f.*Lp.bio;
Mp.die = Lp.con_p.*Lp.bio;

% predation rates (m-2 d-1)
Sf.pred = Sf.die ./ Sf.bio;
Sp.pred = Sp.die ./ Sp.bio;
Mf.pred = Mf.die ./ Mf.bio;
Mp.pred = Mp.die ./ Mp.bio;

% Natural mortality rates
Sf.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Mf.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Lp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Lp.td,M_l);

% Energy available for somatic growth nu
[Sf.nu, Sf.prod] = sub_nu(Sf.I,Sf.bio,Sf.met);
[Sp.nu, Sp.prod] = sub_nu(Sp.I,Sp.bio,Sp.met);
[Mf.nu, Mf.prod] = sub_nu(Mf.I,Mf.bio,Mf.met);
[Mp.nu, Mp.prod] = sub_nu(Mp.I,Mp.bio,Mp.met);
[Lp.nu, Lp.prod] = sub_nu(Lp.I,Lp.bio,Lp.met);

% Maturation (note subscript on Kappa is larvae, juv, adult)
Sf.gamma = sub_gamma(K_l,Z_s,Sf.nu,Sf.die,Sf.bio,Sf.nmort,0,0);
Sp.gamma = sub_gamma(K_l,Z_s,Sp.nu,Sp.die,Sp.bio,Sp.nmort,0,0);
Mf.gamma = sub_gamma(K_a,Z_m,Mf.nu,Mf.die,Mf.bio,Mf.nmort,dfrate,MFsel);
Mp.gamma = sub_gamma(K_j,Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,dfrate,MPsel);
Lp.gamma = sub_gamma(K_a,Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,dfrate,LPsel);

% Egg production (by med and large size classes only)
[Mf.gamma,Mf.nu,Mf.rep,Mf.egg] = sub_rep(Mf.gamma,Mf.nu,K_a,Mf.S(:,DY),Mf.egg);
[Lp.gamma,Lp.nu,Lp.rep,Lp.egg] = sub_rep(Lp.gamma,Lp.nu,K_a,Lp.S(:,DY),Lp.egg);

% Recruitment (from smaller size class)
Sf.rec = sub_rec_larv(Mf.rep,Mf.bio,rfrac);
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,rfrac);
Mf.rec = sub_rec(Sf.gamma,Sf.bio);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);

% Mass balance
Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.egg,Sf.nmort);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.egg,Sp.nmort);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.egg,Mf.nmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.egg,Mp.nmort);

Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.egg,Lp.nmort);

% Fishing by rate
[Mf.bio, Mf.caught, Mf.fmort] = sub_fishing_rate(Mf.bio,dfrate,MFsel);
[Mp.bio, Mp.caught, Mp.fmort] = sub_fishing_rate(Mp.bio,dfrate,MPsel);
[Lp.bio, Lp.caught, Lp.fmort] = sub_fishing_rate(Lp.bio,dfrate,LPsel);

% Advection-Diffusion
% Sf.bio = sub_advec_vel(CGRD,Sf.bio,ENVR.U,ENVR.V,ni,nj,tstep);
% Sp.bio = sub_advec_vel(CGRD,Sp.bio,ENVR.U,ENVR.V,ni,nj,tstep);

% Forward Euler checks for demographics and movement
Sf.bio=sub_check(Sf.bio);
Sp.bio=sub_check(Sp.bio);
Mf.bio=sub_check(Mf.bio);
Mp.bio=sub_check(Mp.bio);
Lp.bio=sub_check(Lp.bio);

end
