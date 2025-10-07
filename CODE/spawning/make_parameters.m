%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function param = make_parameters(param)
    %! Integration parameters
    param.DT = 1.0;       % time step
    param.tstep = 1.0;    % time step in hours for adv-diff

    % define diffusivity
    param.K = 600.0;

    %! Which fishes harvested
    param.MFsel = 1;
    param.LPsel = 1;
    param.LDsel = 1;
    param.Jsel  = 0.1;
    param.MPsel = param.Jsel * param.LPsel;
    param.MDsel = param.Jsel * param.LDsel;

    %! Benthic-pelagic coupling cutoff (depth, m)
    param.PI_be_cutoff = 200;
    % 0:no coupling; 1:demersal coupled only; 2:pelagic & demersal coupled;
    param.pdc = 1;

    %!Individual Mass (g) = geometric mean
    param.M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
    param.M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
    param.M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3
    %logspace(-3,5.0969,7) %gives end points and mid points

    %! Ratio of initial and final body sizes per size-class
    param.Z_s = 0.001/0.5;
    param.Z_m = 0.5/250;
    param.Z_l = 250/125000;

    %%%! Assimilation efficiency lambda (constant across everything)
    param.Lambda = 0.7; %0.579;

    %%%! Kappa rule K as a function of body size
    % K = fraction of energy consumed diverted to somatic growth
    % subscript is larvae, juv, adult)
    param.K_l = 1;
    param.K_j = 1;
    param.K_a = 0.5;

    %%%! Metabolism constants (activity and basal)
    param.amet = 4;       % coeff on met (4, 4.4)
    param.h = 20;         % coeff on Cmax %19.72;
    param.gam = 70;       % coeff on search area (70 baseline) 5.9; 
    param.kc = 0.063;     % coeff on cmax T-dep fn (orig 0.063)
    param.ke = 0.063;     % coeff on enc T-dep fn (orig 0.063)
    param.kt = 0.0855;    % coeff on met T-dep fn (orig 0.063) %0.0855
    param.bpow = 0.175;   % power on metab fn (orig 0.25)
    param.benc = 0.20;    % power on enc fn (orig 0.20)
    param.bcmx = 0.25;    % power on cmax fn (orig 0.25)

    %%%! Transfer efficiency of detritus to benthic prey
    param.bent_eff = 0.075;
    param.CC = 0; % 80

    %%%! Reproductive efficiency
    param.rfrac = 0.01;

    %! Fraction of time spent swimming (from Van Leeuwen)
    param.Tu_s = 1.0;
    param.Tu_m = 1.0; %0.5
    param.Tu_l = 1.0; %0.1

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
    
    param.Sm = 0.25;  %Feeding 2 sizes down, 0.25
    param.D = 0.75;   %Demersal feeding in pelagic reduction
    param.A = 0.5;    %Adult predation reduction %*****

    param.MF_phi_MZ = param.Sm;
    param.MF_phi_LZ = 1.0;
    param.MF_phi_S  = 1.0;

    param.MP_phi_MZ = param.Sm;
    param.MP_phi_LZ = 1.0;
    param.MP_phi_S  = 1.0;

    param.MD_phi_BE = 1.0;

    param.LP_phi_MF = 1.0 * param.A;
    param.LP_phi_MP = 1.0;
    param.LP_phi_MD = 0.0;

    param.LD_phi_MF = param.D * param.A;
    param.LD_phi_MP = param.D;
    param.LD_phi_MD = 1.0;
    param.LD_phi_BE = 1.0;
    
%-----
end
