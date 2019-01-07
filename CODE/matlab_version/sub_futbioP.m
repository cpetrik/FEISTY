%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,ENVR] = sub_futbioP(ID,DY,COBALT,ENVR,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,dfrate,CC)

global DAYS GRD NX
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit 
global bent_eff rfrac Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel LPsel LDsel

%%% COBALT information
ENVR = get_COBALT(COBALT,ID,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);
ENVR.dZm = sub_neg(ENVR.dZm);
ENVR.dZl = sub_neg(ENVR.dZl);

% Pelagic-demersal coupling
%Lp: fraction of time large piscivores spends in pelagic
%Ld: fraction of time large demersals spends in pelagic
if (pdc == 0)
    Lp.td = 1.0;
    Ld.td = 0.0;
elseif (pdc == 1)
    Lp.td = 1.0;
    Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (pdc == 2)
    Lp.td = sub_tdif_pel(ENVR.H,Mf.bio,Mp.bio,Md.bio);
    Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
end

% Metabolism
Sp.met = sub_met(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Mp.met = sub_met(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Lp.met = sub_met(ENVR.Tp,ENVR.Tb,Lp.td,M_l);

% Encounter rates
%sub_enc(Tp,Tb,wgt,prey,td,tprey,pref)
Sp.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sp.td,Sp.td,1);
Mp.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mp.td,Mp.td,MP_phi_MZ);
Mp.enc_zl = sub_enc(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mp.td,Mp.td,MP_phi_LZ);
Mp.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mp.td,Mp.td,MP_phi_S);
Lp.enc_p  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Mp.bio,Lp.td,Lp.td,LP_phi_MP);

% Consumption rates
Sp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sp.td,M_s,Sp.enc_zm);
Mp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zl,Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_d]);
Lp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);

% Offline coupling
%MZ consumption cannot exceed amount lost to higher predation in COBALT runs
[Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
    sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);
%LZ consumption cannot exceed amount lost to higher predation in COBALT runs
[Mf.con_zl,Mp.con_zl,ENVR.fZl] = ...
    sub_offline_zl(Mf.con_zl,Mp.con_zl,Mf.bio,Mp.bio,ENVR.dZl);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sp.I = Sp.con_zm;
Mp.I = Mp.con_zm + Mp.con_zl + Mp.con_p;
Lp.I = Lp.con_p;

% Consumption related to Cmax
Sp.clev = sub_clev(Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Mp.clev = sub_clev(Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Lp.clev = sub_clev(Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,M_l);

% Death rates (g m-2 d-1)
Sp.die = Mp.con_p.*Mp.bio;
Mp.die = Lp.con_p.*Lp.bio;

% predation rates (m-2 d-1)
Sp.pred = Sp.die ./ Sp.bio;
Mp.pred = Mp.die ./ Mp.bio;

% Natural mortality rates
Sp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Mp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Lp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Lp.td,M_l);

% Energy available for somatic growth nu
[Sp.nu, Sp.prod] = sub_nu(Sp.I,Sp.bio,Sp.met);
[Mp.nu, Mp.prod] = sub_nu(Mp.I,Mp.bio,Mp.met);
[Lp.nu, Lp.prod] = sub_nu(Lp.I,Lp.bio,Lp.met);

% Maturation (note subscript on Kappa is larvae, juv, adult)
Sp.gamma = sub_gamma(K_l,Z_s,Sp.nu,Sp.die,Sp.bio,Sp.nmort,0,0);
Mp.gamma = sub_gamma(K_j,Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,0,MPsel);
Lp.gamma = sub_gamma(K_a,Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,dfrate,LPsel);

% Egg production (by med and large size classes only)
[Lp.gamma,Lp.nu,Lp.rep,Lp.egg] = sub_rep(Lp.gamma,Lp.nu,K_a,Lp.S(:,DY),Lp.egg);

% Recruitment (from smaller size class)
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,rfrac);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);

% Mass balance
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.egg,Sp.nmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.egg,Mp.nmort);
Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.egg,Lp.nmort);

% Fishing by rate
[Mp.bio, Mp.caught] = sub_fishing_rate(Mp.bio,dfrate,MPsel);
[Lp.bio, Lp.caught] = sub_fishing_rate(Lp.bio,dfrate,LPsel);

% Forward Euler checks for demographics and movement
Sf.bio=0.0;
Sp.bio=sub_check(Sp.bio);
Sd.bio=0.0;
Mf.bio=0.0;
Mp.bio=sub_check(Mp.bio);
Md.bio=0.0;
Lp.bio=sub_check(Lp.bio);
Ld.bio=0.0;

end
