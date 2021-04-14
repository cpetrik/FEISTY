%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function make_parameters_mzpref(MZpref)
    global D Sm A 
    global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE 
    global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE  
    
    Sm = 0.25;  %Feeding 2 sizes down; 0.25
    D = 0.75;   %Demersal feeding in pelagic reduction
    A = 0.5;    %Adult predation reduction %*****

    MF_phi_MZ = Sm*MZpref;
    MF_phi_LZ = 1.0;
    MF_phi_S  = 1.0;
    
    MP_phi_MZ = Sm*MZpref;
    MP_phi_LZ = 1.0;
    MP_phi_S  = 1.0;
    
    MD_phi_BE = 1.0;
    
    LP_phi_MF = 1.0*A;
    LP_phi_MP = 1.0;
    LP_phi_MD = 0.0;
    
    LD_phi_MF = D*A;
    LD_phi_MP = D;
    LD_phi_MD = 1.0;
    LD_phi_BE = 1.0;
    
%-----
end
