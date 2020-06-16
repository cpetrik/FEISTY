%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sd,Md,Ld,BENT,ENVR] = sub_tbio_benthos_be2(ID,DY,COBALT,ENVR,Sd,Md,Ld,BENT,dfrate)

global NX
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac benc bcmx amet
global Nat_mrt MORT
global MD_phi_BE LD_phi_MD LD_phi_BE
global MDsel LDsel Jsel efn cfn mfn
global kc ke

%%% If biomass < individual fish mass per grid cell, set all rates to zero? %%%

%%% COBALT information
ENVR = get_COBALT(COBALT,ID,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);
ENVR.dZm = sub_neg(ENVR.dZm);
ENVR.dZl = sub_neg(ENVR.dZl);

% Update benthic biomass with new detritus avail at that time step
[BENT.sm,BENT.md,BENT.predS,BENT.predM] = ...
    sub_update_be2_boris(bent_eff,ENVR.det,[BENT.sm,BENT.md],[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);
BENT.sm = sub_check(BENT.sm);
BENT.md = sub_check(BENT.md);

% Pelagic-demersal coupling
%Ld: fraction of time large demersals spends in pelagic
if (pdc == 0)
    Ld.td = zeros(NX,1);
end

% Metabolism
Sd.met = sub_met(ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Md.met = sub_met(ENVR.Tp,ENVR.Tb,Md.td,M_m);
Ld.met = sub_met(ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Encounter rates
%           sub_enc(Tp     ,Tb     ,wgt,prey   ,tpel ,tprey,pref)
Sd.enc_zm = sub_enc(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sd.td,Sd.td,1);

Md.enc_be = sub_enc(ENVR.Tp,ENVR.Tb,M_m,BENT.sm,Md.td,1-Md.td,MD_phi_BE);

Ld.enc_d  = sub_enc(ENVR.Tp,ENVR.Tb,M_l,Md.bio,Ld.td,1-Ld.td,LD_phi_MD);
Ld.enc_be = sub_enc(ENVR.Tp,ENVR.Tb,M_l,BENT.md,Ld.td,1-Ld.td,LD_phi_BE);

% Consumption rates
Sd.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sd.td,M_s,Sd.enc_zm);

Md.con_be = sub_cons(ENVR.Tp,ENVR.Tb,Md.td,M_m,Md.enc_be);

Ld.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_d,Ld.enc_p,Ld.enc_f,Ld.enc_be]);
Ld.con_be = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_be,Ld.enc_f,Ld.enc_p,Ld.enc_d]);

% Offline coupling
%MZ consumption cannot exceed amount lost to higher predation in COBALT runs
[Sd.con_zm,ENVR.fZm] = sub_offline_zm_D(Sd.con_zm,Sd.bio,ENVR.dZm);
%Track fraction of benthic material consumed
[ENVR.fB] = sub_offline_be2([Md.con_be,Ld.con_be],[Md.bio,Ld.bio],BENT.sm,BENT.md);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sd.I = Sd.con_zm;
Md.I = Md.con_be;
Ld.I = Ld.con_f + Ld.con_p + Ld.con_d + Ld.con_be;

% Consumption related to Cmax
Sd.clev = sub_clev(Sd.I,ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Md.clev = sub_clev(Md.I,ENVR.Tp,ENVR.Tb,Md.td,M_m);
Ld.clev = sub_clev(Ld.I,ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Death rates (g m-2 d-1)
%Sd.die = Mp.con_d.*Mp.bio + Mf.con_d.*Mf.bio;
Md.die = Ld.con_d.*Ld.bio;

% predation rates (m-2 d-1)
Sd.pred = Sd.die ./ Sd.bio;
Md.pred = Md.die ./ Md.bio;

% Natural mortality rates
Sd.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Md.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Md.td,M_m);
Ld.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Energy available for somatic growth nu
[Sd.nu, Sd.prod] = sub_nu(Sd.I,Sd.bio,Sd.met);
[Md.nu, Md.prod] = sub_nu(Md.I,Md.bio,Md.met);
[Ld.nu, Ld.prod] = sub_nu(Ld.I,Ld.bio,Ld.met);

% Maturation (note subscript on Kappa is larvae, juv, adult)
Sd.gamma = sub_gamma(K_l,Z_s,Sd.nu,Sd.die,Sd.bio,Sd.nmort,0,0);
Md.gamma = sub_gamma(K_j,Z_m,Md.nu,Md.die,Md.bio,Md.nmort,dfrate,MDsel);
Ld.gamma = sub_gamma(K_a,Z_l,Ld.nu,Ld.die,Ld.bio,Ld.nmort,dfrate,LDsel);

% Egg production (by med and large size classes only)
[Ld.gamma,Ld.nu,Ld.rep,Ld.egg] = sub_rep(Ld.gamma,Ld.nu,K_a,Ld.S(:,DY),Ld.egg);

% Recruitment (from smaller size class)
Sd.rec = sub_rec_larv(Ld.rep,Ld.bio,rfrac);
Md.rec = sub_rec(Sd.gamma,Sd.bio);
Ld.rec = sub_rec(Md.gamma,Md.bio);

% Mass balance
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.egg,Sd.nmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.egg,Md.nmort);
Ld.bio = sub_update_fi(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.egg,Ld.nmort);

% Fishing by rate
[Md.bio, Md.caught, Md.fmort] = sub_fishing_rate(Md.bio,dfrate,MDsel);
[Ld.bio, Ld.caught, Ld.fmort] = sub_fishing_rate(Ld.bio,dfrate,LDsel);

% Forward Euler checks for demographics and movement
Sd.bio=sub_check(Sd.bio);
Md.bio=sub_check(Md.bio);
Ld.bio=sub_check(Ld.bio);

end
