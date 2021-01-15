%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function param = make_parameters_fishvec(param)

%! Integration parameters
param.DT = 1.0;       % time step
param.tstep = 1.0;    % time step in hours for adv-diff

%! Resources:
param.ixR = [1 2 3]; % MZ, LZ, B
param.w(param.ixR) = [2e-06  0.001 0.001];  % lower limit
param.wu(param.ixR) = [0.001 0.5  125];     % upper limit
param.wc(param.ixR) = 10.^((log10(param.w(param.ixR))+log10(param.wu(param.ixR)))./2); %central = geom mean

%! Stages:
param.nstage = 3;       % number of stages predator use 3, 6, 9, etc (prey = 2/3)
param.nsize  = param.nstage + 1; %
param.sfish = 0.001;    % smallest size fish (all fish)
param.lfish = 125000;   % largest size fish (only predator)
param.smat = 0.5;   % weight at maturity forage/meso
param.lmat = 250;   % weight at maturity predators
param.sizes = logspace(log10(param.sfish), log10(param.lfish),param.nsize);
[nb,param.maxsmall] = min(abs(param.sizes-250));
param.stage = param.sizes(end)/param.sizes(end-1);

%! Species:
param.nSpecies = 3;
param.ixFish = 4: 3+(param.nstage*2/3)+(param.nstage)*2; % Index for all fish
% Indices and weight classes small pelagics
param.w = [param.w param.sizes(1:param.maxsmall-1)];
param.ix1(1) = 3+1;
param.ix2(1) = 3+(param.maxsmall-1);
% Indices and weight classes large pelagic:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(2) = param.ix2(1) + 1;
param.ix2(2) = param.ix2(1) + (param.nsize-1);
% Indices and weight classes large demersal:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(3) = param.ix2(2) + 1;
param.ix2(3) = param.ix2(2) + (param.nsize-1);
% Central and upper weights
param.wc(param.ixFish) = param.w(param.ixFish)*sqrt(param.stage); % central sizes
param.wu(param.ixFish) = param.w(param.ixFish)*param.stage; % Upper sizes
%!Individual Mass (g) = geometric mean
% param.M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
% param.M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
% param.M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3

%! Ratio of initial and final body sizes per size-class
param.Z = (param.w./param.wu)';
% param.Z_s = 0.001/0.5;
% param.Z_m = 0.5/250;
% param.Z_l = 250/125000;

%! Benthic-pelagic coupling cutoff (depth, m)
param.PI_be_cutoff = 200;
% 0:no coupling; 1:demersal coupled only; 2:pelagic & demersal coupled;
param.pdc = 1;

%! Define investment in growth vs. repro
[~,param.matstageS] = min(abs(param.sizes-param.smat));
[~,param.matstageL] = min(abs(param.sizes-param.lmat));
param.kappaS = [ones(param.matstageS-1,1)' repmat(0.5,(param.maxsmall-1-(param.matstageS-1)),1)'];
param.kappaL = [ones(param.matstageL-1,1)' repmat(0.5,(param.nsize-1-(param.matstageL-1)),1)'];
param.kappa = [param.kappaS param.kappaL param.kappaL];
%%%! Kappa rule K as a function of body size
% K = fraction of energy consumed diverted to somatic growth
% subscript is larvae, juv, adult)
% param.K_l = 1;
% param.K_j = 1;
% param.K_a = 0.5;

%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;

%! Which fishes harvested
% param.MFsel = 1;
% param.LPsel = 1;
% param.LDsel = 1;
param.Jsel  = 0.1;
param.Asel  = 1;
% param.MPsel = param.Jsel * param.LPsel;
% param.MDsel = param.Jsel * param.LDsel;
param.selS = [0 param.Asel];
param.selL = [0 param.Jsel*param.Asel param.Asel];
param.sel = [param.selS param.selL param.selL];

%! Metabolism constants (activity and basal)
param.Lambda = 0.7;   % Assimilation efficiency lambda (constant across everything)
param.amet = 4;       % coeff on met
param.h = 20;         % coeff on Cmax
param.gam = 70;       % coeff on search area
param.kc = 0.063;     % coeff on cmax T-dep fn (orig 0.063)
param.ke = 0.063;     % coeff on enc T-dep fn (orig 0.063)
param.kt = 0.0855;    % coeff on met T-dep fn (orig 0.063) %0.0855
param.bpow = 0.175;   % power on metab fn (orig 0.25)
param.benc = 0.20;    % power on enc fn (orig 0.20)
param.bcmx = 0.25;    % power on cmax fn (orig 0.25)

%! Transfer efficiency of detritus to benthic prey
param.bent_eff = 0.075;

%! Reproductive efficiency
param.rfrac = 0.01;
param.eRepro = repmat(0.01,param.nSpecies,1)';

%%%! Background mortality
param.Nat_mrt = 0.1/365;
%0=none, 1=constant, 6=const wgt, T-dep, 7=const T, wgt-dep
%2=Hartvig T-dep, 3=mizer T-dep, 4=J&C T-dep, 5=P&W T-dep
param.MORT = 1;

%%%! Diet Preference Phi
% The predator prey mass ratio is assumed 3 orders of mag, i.e. 1000, i.e. one step down
% Because Medium fishes are 2 sizes bigger than Medium zoo, pref = 0.1
% We don't have a pred-prey matrix anymore, we are instead explicit about who eats who:
%-----
%small forage fish eats medium zoo
%small piscivores eats medium zoo
%small detritivore eats medium zoo
%medium forage fish eats medium & large zoo, all small fishes
%medium piscivore eats medium & large zoo, all small fishes
%medium detritivore eats detritus
%large piscivore eats medium forage fish, medium piscivore, medium detritivore
%large detritivore eats detritus, medium forage fish, medium piscivore, medium detrivore

param.Sm = 0.25;  %Feeding 2 sizes down
param.D = 0.75;   %Demersal feeding in pelagic reduction
param.A = 0.75;    %Adult predation reduction %*****

% param.MF_phi_MZ = param.Sm;
% param.MF_phi_LZ = 1.0;
% param.MF_phi_S  = 1.0;
% 
% param.MP_phi_MZ = param.Sm * param.J;
% param.MP_phi_LZ = param.J;
% param.MP_phi_S  = param.J;
% 
% param.MD_phi_BE = 1.0;
% 
% param.LP_phi_MF = 1.0 * param.A;
% param.LP_phi_MP = 1.0;
% param.LP_phi_MD = 0.0;
% 
% param.LD_phi_MF = param.D * param.A;
% param.LD_phi_MP = param.D;
% param.LD_phi_MD = 1.0;
% param.LD_phi_BE = 1.0;

theta = zeros(length(param.wc));

% first stages as juvenile/adult for predators
[~,ixjuv] = min(abs(param.sizes-param.smat));
[~,ixadult] = min(abs(param.sizes-param.lmat));
param.ixS = param.ix1;
param.ixM = param.ix1+1;
param.ixL = param.ix2(2:3);

%! Change specific interactions
% all larvae eat MZ
theta(param.ixS,1) = 1;
% all pelagic medium fish eat MZ at reduced pref
theta(param.ixM(1:2),1) = param.Sm;
% all pelagic medium fish eat LZ
theta(param.ixM(1:2),2) = 1;
% all pelagic medium fish eat Sfish
theta(param.ixM(1:2),param.ixS) = 1;
% pelagic large fish eat Mfish
theta(param.ixL(1),param.ixM(1:2)) = 1;
% demersal large fish eat Mfish at reduced pref
theta(param.ixL(2),param.ixM(1:2)) = param.D;
% demersal large fish eat demersal med
theta(param.ixL(2),param.ixM(3)) = 1;
% benthivory only by demersals
idx_be = param.ixFish(end-1:end); 
theta(idx_be,3) = 1;
% provide benefit to forage fish (predator avoidance)
pred1 = param.ix1(2)+(ixadult-1):param.ix2(2);
pred2 = param.ix1(3)+(ixadult-1):param.ix2(3);
prey1 = param.ix1(1)+(ixjuv-1):param.ix2(1);
idx_predat = [pred1 pred2];
idx_prey= prey1;
theta(idx_predat,idx_prey) = theta(idx_predat,idx_prey)*param.A;

param.theta = theta;

%-----
end
