%%%%!! RUN Climatol FOR ALL LOCATIONS
clear all
close all

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj frate kc ke

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
Pdrpbx = '/Users/cpetrik/Dropbox/';
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat');
NX = length(GRD.Z);
ID = 1:NX;

%%%%%%%%%%%%%%% Initialize Model Variables
%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

paramt = cell(39,1);
simt = cell(39,1);
% PARAMETER SENSITIVITY TEST
for j = 1:39
    
    %! Make core parameters/constants (global)
    make_parameters_mid()

    %! Change individual parameters
    if j ==1
        h = 5;
        ptext = 'h5';
    elseif j==2
        h = 500;
        ptext = 'h500';
    elseif j==3
        gam = 5;
        ptext = 'gam5';
    elseif j==4
        gam = 500;
        ptext = 'gam500';
    elseif j==5
        amet = 2;
        ptext = 'amet2';
    elseif j==6
        amet = 8;
        ptext = 'amet8';
    elseif j==7
        Lambda = 0.56;
        ptext = 'lam056';
    elseif j==8
        Lambda = 0.84;
        ptext = 'lam084';
    elseif j==9
        bcmx = 0.1;
        ptext = 'bc100';
    elseif j==10
        bcmx = 0.32;
        ptext = 'bc320';
    elseif j==11
        benc = 0.1;
        ptext = 'be100';
    elseif j==12
        benc = 0.32;
        ptext = 'be320';
    elseif j==13
        bpow = 0.1;
        ptext = 'bm100';
    elseif j==14
        bpow = 0.32;
        ptext = 'bm320';
    elseif j==15
        bent_eff = 0.0375;
        ptext = 'BE0375';
    elseif j==16
        bent_eff = 0.15;
        ptext = 'BE15';
    elseif j==17
        rfrac = 0.001;
        ptext = 'RE0001';
    elseif j==18
        rfrac = 0.1;
        ptext = 'RE01';
    elseif j==19
        frate = 0.15;
        ptext = 'fish015';
    elseif j==20
        frate = 0.6;
        ptext = 'fish06';
    elseif j==21
        kc = 0.0302;
        ptext = 'kc0302';
    elseif j==22
        kc = 0.1208;
        ptext = 'kc1208';
    elseif j==23
        ke = 0.0302;
        ptext = 'ke0302';
    elseif j==24
        ke = 0.1208;
        ptext = 'ke1208';
    elseif j==25
        kt = 0.0302;
        ptext = 'kt0302';
    elseif j==26
        kt = 0.1208;
        ptext = 'kt1208';
    elseif j==27
        K_a = 0.25;
        ptext = 'kap25';
    elseif j==28
        K_a = 0.75;
        ptext = 'kap75';
    elseif j==29
        Nat_mrt = 0.05/365;
        ptext = 'unat05';
    elseif j==30
        Nat_mrt = 0.2/365;
        ptext = 'unat20';
    elseif j==31
        A = 0.25;
        ptext = 'A025';
    elseif j==32
        A = 1.0;
        ptext = 'A100';
    elseif j==33
        Sm = 0.125;
        ptext = 'Sm0125';
    elseif j==34
        Sm = 0.5;
        ptext = 'Sm050';
    elseif j==35
        J = 0.25;
        ptext = 'J025';
    elseif j==36
        J = 1.0;
        ptext = 'J100';
    elseif j==37
        D = 0.25;
        ptext = 'D025';
    elseif j==38
        D = 1.0;
        ptext = 'D100';
    else
        ptext = 'base';
    end
      
    %! Set globals that are function of other parameters
    dfrate = frate/365.0;
    MF_phi_MZ = Sm;
    MP_phi_MZ = Sm*J;
    MP_phi_LZ = J;
    MP_phi_S  = J;
    LP_phi_MF = 1.0*A;
    LD_phi_MF = D*A;
    LD_phi_MP = D;
    
    %! Create a directory for output
    fname = sub_fname_param_v2(ptext);
    oname = sub_fname_nodir(frate);
        
    paramt{j} = ptext;
    simt{j} = oname;
    
end

cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];
save([nfile 'Climatol_All_fish03_means_param_sens_v2.mat'],'paramt','simt',...
    '-append')

