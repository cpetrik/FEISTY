%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function make_params_benthos_be2()
    global DAYS GRD NX ID
    global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l
    global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
    global bent_eff rfrac benc bcmx amet
    global Nat_mrt MORT
    global MD_phi_BE LD_phi_MD LD_phi_BE
    global MDsel LDsel Jsel efn cfn mfn
    global kc ke 

    %! Integration parameters
    DT = 1.0;       % time step
    tstep = 1.0;    % time step in hours for adv-diff
    
    %! Which fishes harvested
    LDsel = 1;
    Jsel  = 0.1;
    MDsel = Jsel * LDsel;
    
    %! Benthic-pelagic coupling cutoff (depth, m)
    PI_be_cutoff = 200;
    % 0:no coupling; 1:demersal coupled only; 2:pelagic & demersal coupled;
    pdc = 0;

    %!Individual Mass (g) = geometric mean
    M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
    M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
    M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3
    %logspace(-3,5.0969,7) %gives end points and mid points

    %! Body lengths (mm)
    % Convert from mm to cm and use their const coeff = 0.01g/cm3
    L_s = 10.0 * (M_s/0.01)^(1/3); % small   4.6-36.8 mm, 13.1 mm
    L_m = 10.0 * (M_m/0.01)^(1/3); % medium 36.8-292 mm, 10.4 cm
    L_l = 10.0 * (M_l/0.01)^(1/3); % large  0.292-2.32 m, 0.82 m

    %! Ratio of initial and final body sizes per size-class
    Z_s = 0.001/0.5;
    Z_m = 0.5/250;
    Z_l = 250/125000;

    %%%! Assimilation efficiency lambda (constant across everything)
    Lambda = 0.7;

    %%%! Kappa rule K as a function of body size
    % K = fraction of energy consumed diverted to somatic growth
    % subscript is larvae, juv, adult)
    K_l = 1;
    K_j = 1;
    K_a = 0.5;

    %%%! Metabolism constants (activity and basal)
    amet = 4;       % coeff on met
    h = 20;         % coeff on Cmax
    gam = 70;       % coeff on search area
    kc = 0.063;     % coeff on cmax T-dep fn (orig 0.063)
    ke = 0.063;     % coeff on enc T-dep fn (orig 0.063)
    kt = 0.0855;    % coeff on met T-dep fn (orig 0.063) %0.0855
    bpow = 0.175;   % power on metab fn (orig 0.25)
    benc = 0.20;    % power on enc fn (orig 0.20)
    bcmx = 0.25;    % power on cmax fn (orig 0.25)
    
    %%%! Benthic prey
    bent_eff = 0.75;   % Transfer of detritus to benthos (not bacteria)
    
    %%%! Reproductive efficiency
    rfrac = 0.01;

    %%%! Background mortality
    Nat_mrt = 0.1/365; 
    %0=none, 1=constant, 6=const wgt, T-dep, 7=const T, wgt-dep
    %2=Hartvig T-dep, 3=mizer T-dep, 4=J&C T-dep, 5=P&W T-dep
    MORT = 1;   

    %%%! Diet Preference Phi
    % The predator prey mass ratio is assumed 3 orders of mag, i.e. 1000, i.e. one step down
    % Because Medium fishes are 2 sizes bigger than Medium zoo, pref = 0.1
    % We don't have a pred-prey matrix anymore, we are instead explicit about who eats who:
    %-----
    %small detritivore eats medium zoo
    %medium detritivore eats detritus
    %large detritivore eats detritus, medium detrivore
    
    MD_phi_BE = 1.0;
    LD_phi_MD = 1.0;
    LD_phi_BE = 1.0;
    
%-----
end
