%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [fish,BENT,ENVR] = sub_futbio_fishvec(ID,DY,COBALT,GRD,fish,BENT,param)

%%% If biomass < individual fish mass per grid cell, set all rates to zero? %%%

%%% COBALT information
ENVR = get_COBALT(COBALT,GRD,ID,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.Zl  = sub_neg(ENVR.Zl);
ENVR.dZm = sub_neg(ENVR.dZm);
ENVR.dZl = sub_neg(ENVR.dZl);
fish.bio(:,1) = ENVR.Zm;
fish.bio(:,2) = ENVR.Zl;

% Update benthic biomass with new detritus avail at that time step
[BENT.mass,BENT.pred] = sub_update_be_fishvec(BENT.mass,param,ENVR.det,fish.con,fish.bio);
BENT.mass = sub_check(BENT.mass);
fish.bio(:,3) = BENT.mass;

% Pelagic-demersal coupling
%Lp: fraction of time large piscivores spends in pelagic
%Ld: fraction of time large demersals spends in pelagic
if (param.pdc == 0)
    fish.td(:,8)  = ones(param.NX,1);
    fish.td(:,11) = zeros(param.NX,1);
elseif (param.pdc == 1)
    fish.td(:,8)  = ones(param.NX,1);
    fish.td(:,11) = sub_tdif_dem_fishvec(ENVR.H,param,fish.bio,BENT.mass);
elseif (param.pdc == 2)
    fish.td(:,8)  = sub_tdif_pel_fishvec(ENVR.H,param,fish.bio);
    fish.td(:,11) = sub_tdif_dem_fishvec(ENVR.H,param,fish.bio,BENT.mass);
else
    fish.td(:,8) = ones(param.NX,1);
    fish.td(:,11) = ones(param.NX,1);
end

% Calc temp
Tp = repmat(ENVR.Tp,1,length(param.wc));
Tb = repmat(ENVR.Tb,1,length(param.wc));
temp = (Tp .* fish.td) + (Tb .* (1.0-fish.td));

% Metabolism
fish.met = (exp(param.kt * (temp-10.0)) .* param.amet .* param.wc.^(-param.bpow)) ./365.0;

% Encounter rates
%Search rate
A = (exp(param.ke * (temp-10.0)) .* param.gam .* param.wc.^(-param.benc)) ./365.0;
%dim(enc) = [NX,pred,prey]
fish.enc = sub_enc_fishvec(param,A,fish);

% Consumption rates
%Cmax rate
cmax = (exp(param.kc * (temp-10.0)) .* param.h .* param.wc.^(-param.bcmx)) ./365.0;
%Type II consumption
fish.con = sub_cons_fishvec(cmax,fish);

% Offline coupling
%MZ consumption cannot exceed amount lost to higher predation in COBALT runs
[fish,ENVR] = sub_off_zm_fishvec(fish,ENVR);
%LZ consumption cannot exceed amount lost to higher predation in COBALT runs
[fish,ENVR] = sub_off_zl_fishvec(fish,ENVR);
%Track fraction of benthic material consumed
ENVR = sub_off_be_fishvec(fish,ENVR);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
fish.I = sum(fish.con,3);

% Consumption related to Cmax ("feeding level")
fish.clev = fish.I ./ cmax;

% Death rates (g m-2 d-1)
fish.die = sub_pred_fishvec(param,fish);

% biomass-specific predation rates (m-2 d-1)
fish.pred = fish.die ./ fish.bio;

% Natural mortality rates
fish.nmort = sub_nmort_fishvec(param,fish,temp);

% Energy available for somatic growth nu
fish.nu = (fish.I .* param.Lambda) - fish.met;
fish.prod = fish.nu .* fish.bio;

% Maturation (note subscript on Kappa is larvae, juv, adult)
gamma = sub_gamma_fishvec(param,fish);
fish.gamma(:,param.ixFish) = gamma;

% Egg production (by adults only)
[gamma,rep] = sub_rep_fishvec(param,fish);
fish.gamma(:,param.ixFish) = gamma;
fish.rep(:,param.ixFish) = rep;

% Recruitment (biomass, not rate)
%larvae are "recruited" through reproduction of adults
fish.rec(:,param.ix1) = param.rfrac .* fish.rep(:,param.ix2) .* fish.bio(:,param.ix2);
%recruitment from smaller size class
ix = setdiff(param.ixFish,param.ix1);
fish.rec(:,ix) = fish.gamma(:,ix-1) .* fish.bio(:,ix-1);

% Mass balance
% all inputs are in g g-1 d-1, except rec & die are g d-1
db = fish.rec + ((fish.nu - fish.rep - fish.gamma - fish.nmort) .* fish.bio) - fish.die;
fish.bio(:,param.ixFish) =  fish.bio(:,param.ixFish) + db(:,param.ixFish);

% Fishing by rate
fish.caught(:,param.ixFish) = fish.bio(:,param.ixFish) .* param.sel .* param.dfrate;
fish.fmort = fish.caught ./ fish.bio;
fish.bio = fish.bio - fish.caught;

% Forward Euler checks for demographics and movement
fish.bio(fish.bio<0) = eps();

end
