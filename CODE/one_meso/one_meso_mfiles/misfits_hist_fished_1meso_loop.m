% Loop over output of assim & amet values

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

%%
param.frate = 0.3;
param.dfrate = param.frate/365.0;

lam = 0.55:0.025:0.675;
aM = 4.25:0.25:6;

nfiles = cell(length(lam)*length(aM),1);
sfiles = cell(length(lam)*length(aM),1);

i=0;
for L = 1:length(lam)
    for m = 1:length(aM)
        
        i=i+1;
        
        %! Make core parameters/constants (global)
        param = make_parameters_1meso(param); % make core parameters/constants
        
        param.Lambda = lam(L);
        param.amet = aM(m);
        
        %! Create a directory for output
        [fname,simname] = sub_fname_hist_1meso(param);
        nfiles{i} = fname;
        sfiles{i} = simname;
        
    end
end

%%
mis_all = NaN*ones(length(nfiles),48111,8);
mis_sau = NaN*ones(length(nfiles),45,5);

%%
for j=1:length(nfiles)
    load([nfiles{j},'_Means.mat']);
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
    
    %% lme
    load([dpath 'LME_Hist_1meso_All_fish03_.mat'])
    % didn't save yield, try 0.3 * biomass
%     Frate = 0.3/365;
%     mp_my50 = 0.1*Frate * mp_mean50;
%     md_my50 = 0.1*Frate * md_mean50;
%     mf_my50 = Frate * mf_mean50;
%     lp_my50 = Frate * lp_mean50;
%     ld_my50 = Frate * ld_mean50;
%     
%     [lme_mcatch,lme_area] = lme_hist_fished_1meso_loop(mf_my50,...
%     mp_my50,md_my50,lp_my50,ld_my50,dpath);

    %% saup comp
    [r,rmse,ss,mis] = hist_1meso_lme_saup_corr_loop(lme_mcatch,lme_area,ppath);
    
    %% diff from 1 meso w/overcon
    [misf] = hist_1meso_ratios_misfit_loop(sF,sP,sD,sB,sM,sL);
    
    %% save as matrices
    mis_all(j,:,:) = misf;
    mis_sau(j,:,:) = mis;
    
end

%%
save([dp 'bio_rates/Lambda_amet_search_1meso_Dc_enc70-b200_m-b175-k086_c20-b250_D080_A050_nmort1_BE10_noCC_RE00100.mat'],...
    'mis_all','mis_sau','-append')
