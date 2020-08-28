% Loop over output of aE, BE, and A values

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

%%
param.frate = 0.3;
param.dfrate = param.frate/365.0;

aE = 5:8;
BE = [0.05,0.075];
A = [0.5,0.75];

nfiles = cell(4*length(aE),1);
sfiles = cell(4*length(aE),1);
pset = NaN*ones(4*length(aE),4);
%%
i=0;
for x = 1:length(aE)
    for y = 1:length(BE)
        for z = 1:length(A)
            
            i=i+1;
            
            %! Make core parameters/constants (global)
            param = make_parameters_1meso(param); % make core parameters/constants
            
            param.gam = aE(x);
            param.bent_eff = BE(y);
            param.A = A(z);
            
            pset(i,1) = param.h;
            pset(i,2) = param.gam;
            pset(i,3) = param.bent_eff;
            pset(i,4) = param.A;
            
            %! Create a directory for output
            [fname,simname] = sub_fname_hist_1meso(param);
            nfiles{i} = fname;
            sfiles{i} = simname;
            
        end
    end
end

%% manually add a few more combinations
%enc = [7,7]; a = [0.5,0.5];
be = [0.075,0.1];
for n=1:2
    param = make_parameters_1meso(param); % make core parameters/constants
    
    param.gam = 7;
    param.h = 19;
    param.bent_eff = be(n);
    param.A = 0.5;
    
    pset(i+n,1) = param.h;
    pset(i+n,2) = param.gam;
    pset(i+n,3) = param.bent_eff;
    pset(i+n,4) = param.A;
    
    %! Create a directory for output
    [fname,simname] = sub_fname_hist_1meso(param);
    nfiles{i+n} = fname;
    sfiles{i+n} = simname;
end

%%
mis_all = NaN*ones(length(nfiles),48111,9);
mis_sau = NaN*ones(length(nfiles),45,5);

%%
nr = [2,6,10,14];
ran = [1,3:5,7:9,11:13,15:18];
for k=1:length(ran)%1:length(nfiles)
    j = ran(k);
    
    ftex = [nfiles{j}(1:114) 'Means_Historic_1meso_All_fish03_' sfiles{j} '.mat'];
    load(ftex);
    
    dpath = [dp sfiles{j} '/'];
    ppath = [pp sfiles{j} '/'];
    if (~isfolder(ppath))
        mkdir(ppath)
    end
    
    %% space means
    sF = sf_mean50+mf_mean50;
    sP = sp_mean50+mp_mean50+lp_mean50;
    sD = sd_mean50+md_mean50+ld_mean50;
    sM = mf_mean50+mp_mean50+md_mean50;
    sL = lp_mean50+ld_mean50;
    sB = b_mean50;
    sZ = mz_mfrac5;
    
    %% lme
    load([dpath 'LME_Hist_1meso_All_fish03_' sfiles{j} '.mat'])
    
    %% saup comp
    [r,rmse,ss,mis] = hist_1meso_lme_saup_corr_loop(lme_mcatch,lme_area,ppath);
    
    %% diff from 1 meso w/overcon
    [misf] = hist_1meso_ratios_zoo_misfit_loop(sF,sP,sD,sB,sM,sL,sZ);
    
    %% save as matrices
    mis_all(j,:,:) = misf;
    mis_sau(j,:,:) = mis;
    
end

%%
save([dp 'bio_rates/aEnc_BE_A_search_1meso_Dc_enc-b200_m4-b175-k086_c20-b250_D080_nmort1_noCC_RE00100.mat'],...
    'mis_all','mis_sau','nfiles','sfiles','pset')%,'-append')
