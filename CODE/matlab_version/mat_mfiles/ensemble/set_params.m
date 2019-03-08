%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function set_params(X)
    global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl 
    global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
    global bent_eff rfrac D J Sm A benc bcmx amet 
    global Tu_s Tu_m Tu_l Nat_mrt MORT
    global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE 
    global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE  
    global MFsel MPsel MDsel LPsel LDsel Jsel
    global tstep K frate dfrate kc ke
    
    %! Set fishing rate
    frate = X(1); %Fish(F);
    dfrate = frate/365.0;
    
    %%%! Assimilation efficiency lambda (constant across everything)
    Lambda = X(2);

    %%%! Kappa rule K as a function of body size
    % K = fraction of energy consumed diverted to somatic growth
    % subscript is larvae, juv, adult)
    K_a = X(3);

    %%%! Metabolism constants (activity and basal)
    amet = X(4);       % coeff on met
    h    = X(5);         % coeff on Cmax
    gam  = X(6);       % coeff on search area
    kc = X(7);     % coeff on cmax T-dep fn 
    ke = X(8);     % coeff on enc T-dep fn 
    kt = X(9);    % coeff on met T-dep fn (orig 0.063)
    bpow = X(10);   % power on metab fn (orig 0.25)
    benc = X(11);    % power on enc fn (orig 0.20)
    bcmx = X(12);    % power on cmax fn (orig 0.25)
    
    %%%! Transfer efficiency of detritus to benthic prey
    bent_eff = X(13);
    
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
    
    J = X(14);    %Juvenile feeding reduction
    A = X(15);    %Adult predation reduction %*****

    
%-----
end
