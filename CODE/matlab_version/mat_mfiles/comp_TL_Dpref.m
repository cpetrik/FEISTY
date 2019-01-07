% Compare TEs and TLs

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';
coast = [3;4;6;8;10];
coastal = spots(coast)';

% POEM file info
tharv = 'Harvest all fish 0.3 yr^-^1';
sname = 'Clim_locs_All_fish03';
harv = 'All_fish03';

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%% Zoop and det
cpath = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
gpath='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([gpath 'clim_det_biom_Dmeans_Ytot.mat'])
load([cpath 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%ESM2.6 in mg C m-2 or mg C m-2 d-1
%from mg C m-2 to g(WW) m-2
% 1e-3 g C in 1 mg C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)

mmz_mean = mz_mean_clim(ID) * 1e-3 * 9.0;
mlz_mean = lz_mean_clim(ID) * 1e-3 * 9.0;
mmz_loss = mzloss_mean_clim(ID) * 1e-3 * 9.0;
mlz_loss = lzloss_mean_clim(ID) * 1e-3 * 9.0;

tmz_mean = mz_tot_clim(ID) * 1e-3 * 9.0;
tlz_mean = lz_tot_clim(ID) * 1e-3 * 9.0;
tmz_loss = mzl_tot_clim(ID) * 1e-3 * 9.0;
tlz_loss = lzl_tot_clim(ID) * 1e-3 * 9.0;

mdet = det_mean_clim(ID)* 1e-3 * 9.0;
tdet = det_tot_clim(ID)* 1e-3 * 9.0;

mmz_locs = mmz_mean(ids);
mlz_locs = mlz_mean(ids);
tmz_locs = tmz_mean(ids);
tlz_locs = tlz_mean(ids);
mmzl_locs = mmz_loss(ids);
mlzl_locs = mlz_loss(ids);
tmzl_locs = tmz_loss(ids);
tlzl_locs = tlz_loss(ids);
mdet_locs = mdet(ids);
tdet_locs = tdet(ids);

%%
dees = 0.1:0.1:1;
TL = NaN*ones(length(dees),length(spots));
TEc = TL;
TEp1 = TL;
TEp2 = TL;
TEe1 = TL;
TEe2 = TL;
for e=1:length(dees)
    D = dees(e);
    td = num2str(1000+int64(100*D));
    close all
    dp = ['Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D',td(2:end),...
        '_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
    dpath = [datap char(dp) '/'];
    fpath = [figp char(dp) '/'];
    cfile = char(dp);
    load([dpath sname '.mat'])
    
    %%
    all_mean=NaN*ones(3,4,length(spots));
    z = NaN*ones(length(spots),3);
    
    SP = S_Sml_p;
    SF = S_Sml_f;
    SD = S_Sml_d;
    MP = S_Med_p;
    MF = S_Med_f;
    MD = S_Med_d;
    LP = S_Lrg_p;
    LD = S_Lrg_d;
    CO = S_Cobalt;
    
    t=1:size(SP,1);
    lyr=t((end-12+1):end);
    
    %%
    prey = NaN*ones(8,6,length(spots));
    TL1 = prey;
    %fixed TL of prey items
    %should benthic inverts be >2?
    tl(1,:) = [2,2.5,3,3,3,2];
    tl(2,:) = [2,2.5,3,3,3,2];
    tl(3,:) = [2,2.5,3,3,3,2];
    tl(4,:) = [2,2.5,3,3,3,2];
    tl(5,:) = [2,2.5,3,3,3,2];
    tl(6,:) = [2,2.5,3,3,3,2];
    %initial guess for medium TL
    tl(7,:) = [2,2.5,3.5,3.5,3,2];
    tl(8,:) = [2,2.5,3.5,3.5,3,2];
    
    tl_mat = repmat(tl,1,1,16);
    
    %% Final mean biomass in each size
    SP_mean=squeeze(mean(SP(lyr,1,:)))';
    
    %% Growth rate (nu - energy for biomass production)
    SP_mgr=squeeze(nanmean(SP(lyr,15,:)))';
    
    %% Consump per biomass (I) by type
    conF(:,1)=squeeze(nanmean(SF(lyr,8,:)))+squeeze(nanmean(MF(lyr,8,:)));
    conF(:,2)=squeeze(nanmean(SP(lyr,8,:)))+squeeze(nanmean(MP(lyr,8,:)))+squeeze(nanmean(LP(lyr,8,:)));
    conF(:,3)=squeeze(nanmean(SD(lyr,8,:)))+squeeze(nanmean(MD(lyr,8,:)))+squeeze(nanmean(LD(lyr,8,:)));
    conP(:,1)=squeeze(nanmean(SF(lyr,9,:)))+squeeze(nanmean(MF(lyr,9,:)));
    conP(:,2)=squeeze(nanmean(SP(lyr,9,:)))+squeeze(nanmean(MP(lyr,9,:)))+squeeze(nanmean(LP(lyr,9,:)));
    conP(:,3)=squeeze(nanmean(SD(lyr,9,:)))+squeeze(nanmean(MD(lyr,9,:)))+squeeze(nanmean(LD(lyr,9,:)));
    conD(:,1)=squeeze(nanmean(SF(lyr,10,:)))+squeeze(nanmean(MF(lyr,10,:)));
    conD(:,2)=squeeze(nanmean(SP(lyr,10,:)))+squeeze(nanmean(MP(lyr,10,:)))+squeeze(nanmean(LP(lyr,10,:)));
    conD(:,3)=squeeze(nanmean(SD(lyr,10,:)))+squeeze(nanmean(MD(lyr,10,:)))+squeeze(nanmean(LD(lyr,10,:)));
    conZm(:,1)=squeeze(nanmean(SF(lyr,11,:)))+squeeze(nanmean(MF(lyr,11,:)));
    conZm(:,2)=squeeze(nanmean(SP(lyr,11,:)))+squeeze(nanmean(MP(lyr,11,:)))+squeeze(nanmean(LP(lyr,11,:)));
    conZm(:,3)=squeeze(nanmean(SD(lyr,11,:)))+squeeze(nanmean(MD(lyr,11,:)))+squeeze(nanmean(LD(lyr,11,:)));
    conZl(:,1)=squeeze(nanmean(SF(lyr,12,:)))+squeeze(nanmean(MF(lyr,12,:)));
    conZl(:,2)=squeeze(nanmean(SP(lyr,12,:)))+squeeze(nanmean(MP(lyr,12,:)))+squeeze(nanmean(LP(lyr,12,:)));
    conZl(:,3)=squeeze(nanmean(SD(lyr,12,:)))+squeeze(nanmean(MD(lyr,12,:)))+squeeze(nanmean(LD(lyr,12,:)));
    conB(:,1)=squeeze(nanmean(SF(lyr,13,:)))+squeeze(nanmean(MF(lyr,13,:)));
    conB(:,2)=squeeze(nanmean(SP(lyr,13,:)))+squeeze(nanmean(MP(lyr,13,:)))+squeeze(nanmean(LP(lyr,13,:)));
    conB(:,3)=squeeze(nanmean(SD(lyr,13,:)))+squeeze(nanmean(MD(lyr,13,:)))+squeeze(nanmean(LD(lyr,13,:)));
    
    %% TL calc from prey items
    prey(1,1,:)=squeeze(nanmean(SF(lyr,11,:)))';
    prey(1,2,:)=squeeze(nanmean(SF(lyr,12,:)))';
    prey(1,3,:)=squeeze(nanmean(SF(lyr,8,:)))';
    prey(1,4,:)=squeeze(nanmean(SF(lyr,9,:)))';
    prey(1,5,:)=squeeze(nanmean(SF(lyr,10,:)))';
    prey(1,6,:)=squeeze(nanmean(SF(lyr,13,:)))';
    
    prey(2,1,:)=squeeze(nanmean(SP(lyr,11,:)))';
    prey(2,2,:)=squeeze(nanmean(SP(lyr,12,:)))';
    prey(2,3,:)=squeeze(nanmean(SP(lyr,8,:)))';
    prey(2,4,:)=squeeze(nanmean(SP(lyr,9,:)))';
    prey(2,5,:)=squeeze(nanmean(SP(lyr,10,:)))';
    prey(2,6,:)=squeeze(nanmean(SP(lyr,13,:)))';
    
    prey(3,1,:)=squeeze(nanmean(SD(lyr,11,:)))';
    prey(3,2,:)=squeeze(nanmean(SD(lyr,12,:)))';
    prey(3,3,:)=squeeze(nanmean(SD(lyr,8,:)))';
    prey(3,4,:)=squeeze(nanmean(SD(lyr,9,:)))';
    prey(3,5,:)=squeeze(nanmean(SD(lyr,10,:)))';
    prey(3,6,:)=squeeze(nanmean(SD(lyr,13,:)))';
    
    prey(4,1,:)=squeeze(nanmean(MF(lyr,11,:)))';
    prey(4,2,:)=squeeze(nanmean(MF(lyr,12,:)))';
    prey(4,3,:)=squeeze(nanmean(MF(lyr,8,:)))';
    prey(4,4,:)=squeeze(nanmean(MF(lyr,9,:)))';
    prey(4,5,:)=squeeze(nanmean(MF(lyr,10,:)))';
    prey(4,6,:)=squeeze(nanmean(MF(lyr,13,:)))';
    
    prey(5,1,:)=squeeze(nanmean(MP(lyr,11,:)))';
    prey(5,2,:)=squeeze(nanmean(MP(lyr,12,:)))';
    prey(5,3,:)=squeeze(nanmean(MP(lyr,8,:)))';
    prey(5,4,:)=squeeze(nanmean(MP(lyr,9,:)))';
    prey(5,5,:)=squeeze(nanmean(MP(lyr,10,:)))';
    prey(5,6,:)=squeeze(nanmean(MP(lyr,13,:)))';
    
    prey(6,1,:)=squeeze(nanmean(MD(lyr,11,:)))';
    prey(6,2,:)=squeeze(nanmean(MD(lyr,12,:)))';
    prey(6,3,:)=squeeze(nanmean(MD(lyr,8,:)))';
    prey(6,4,:)=squeeze(nanmean(MD(lyr,9,:)))';
    prey(6,5,:)=squeeze(nanmean(MD(lyr,10,:)))';
    prey(6,6,:)=squeeze(nanmean(MD(lyr,13,:)))';
    
    prey(7,1,:)=squeeze(nanmean(LP(lyr,11,:)))';
    prey(7,2,:)=squeeze(nanmean(LP(lyr,12,:)))';
    prey(7,3,:)=squeeze(nanmean(LP(lyr,8,:)))';
    prey(7,4,:)=squeeze(nanmean(LP(lyr,9,:)))';
    prey(7,5,:)=squeeze(nanmean(LP(lyr,10,:)))';
    prey(7,6,:)=squeeze(nanmean(LP(lyr,13,:)))';
    
    prey(8,1,:)=squeeze(nanmean(LD(lyr,11,:)))';
    prey(8,2,:)=squeeze(nanmean(LD(lyr,12,:)))';
    prey(8,3,:)=squeeze(nanmean(LD(lyr,8,:)))';
    prey(8,4,:)=squeeze(nanmean(LD(lyr,9,:)))';
    prey(8,5,:)=squeeze(nanmean(LD(lyr,10,:)))';
    prey(8,6,:)=squeeze(nanmean(LD(lyr,13,:)))';
    
    %% Trophic level calc
    wgt_prey = prey ./ repmat(sum(prey,2),1,6);
    TL1(1:6,:,:) = tl_mat(1:6,:,:) .* wgt_prey(1:6,:,:);
    TLmed = sum(TL1(4:6,:,:),2) + 1;
    tl_mat(7,3:5,:) = TLmed;
    tl_mat(8,3:5,:) = TLmed;
    TL1(7:8,:,:) = tl_mat(7:8,:,:) .* wgt_prey(7:8,:,:);
    TLprey = squeeze(sum(TL1,2)) + 1;
    
    %% Consump per biomass (I)
    SP_con=squeeze(nanmean(SP(lyr,14,:)))';
    SF_con=squeeze(nanmean(SF(lyr,14,:)))';
    SD_con=squeeze(nanmean(SD(lyr,14,:)))';
    MP_con=squeeze(nanmean(MP(lyr,14,:)))';
    MF_con=squeeze(nanmean(MF(lyr,14,:)))';
    MD_con=squeeze(nanmean(MD(lyr,14,:)))';
    LP_con=squeeze(nanmean(LP(lyr,14,:)))';
    LD_con=squeeze(nanmean(LD(lyr,14,:)))';
    
    Pcon=[SP_con;MP_con;LP_con];
    Fcon=[SF_con;MF_con];
    Dcon=[SD_con;MD_con;LD_con];
    
    %% Consumption efficiency
    % Biomass consumed (g/day) (I * biom)
    S_con=(SP_con+SF_con+SD_con);
    M_con=(MP_con+MF_con+MD_con);
    L_con=(LP_con+LD_con);
    
    Con = [S_con;M_con;L_con];
    
    TEcon       = (M_con./S_con);
    TEcon(2,:)  = (L_con./M_con);
    
    %% Production (= nu * biom)
    SP_prod=squeeze(nanmean(SP(lyr,21,:)))';
    SF_prod=squeeze(nanmean(SF(lyr,21,:)))';
    SD_prod=squeeze(nanmean(SD(lyr,21,:)))';
    MP_prod=squeeze(nanmean(MP(lyr,21,:)))';
    MF_prod=squeeze(nanmean(MF(lyr,21,:)))';
    MD_prod=squeeze(nanmean(MD(lyr,21,:)))';
    LP_prod=squeeze(nanmean(LP(lyr,21,:)))';
    LD_prod=squeeze(nanmean(LD(lyr,21,:)))';
    B_prod =squeeze(nanmean(CO(lyr,1,:)))';
    
    Pprod=[SP_prod;MP_prod;LP_prod];
    Fprod=[SF_prod;MF_prod];
    Dprod=[SD_prod;MD_prod;LD_prod];
    
    %% Production (= nu * biom) efficiency
    %Try using losses vs. standing biomass
    MZ_prod = mmz_locs';
    LZ_prod = mlz_locs';
    MZl_prod = mmzl_locs';
    LZl_prod = mlzl_locs';
    
    S_prod=(SP_prod+SF_prod+SD_prod);
    M_prod=(MP_prod+MF_prod+MD_prod);
    L_prod=(LP_prod+LD_prod);
    
    Prod = [MZl_prod;LZl_prod;B_prod;S_prod;M_prod;L_prod];
    
    TEprod_grp      = S_prod./MZl_prod;
    TEprod_grp(2,:) = M_prod./S_prod;
    TEprod_grp(3,:) = L_prod./M_prod;
    
    M_prod1=(MP_prod+MF_prod);
    M_prod2=(MD_prod);
    L_prod1=(LP_prod);
    L_prod2=(LD_prod);
    
    TEprod      = S_prod./MZl_prod;
    TEprod(2,:) = M_prod1./(S_prod+LZl_prod);
    TEprod(3,:) = M_prod2./(B_prod);
    TEprod(4,:) = L_prod1./M_prod1;
    TEprod(5,:) = L_prod2./(M_prod2+B_prod);
    
    %% Effective TE
    TEeff1      = M_prod./(B_prod + MZ_prod + LZ_prod);
    TEeff1(2,:) = L_prod./(B_prod + MZ_prod + LZ_prod);
    
    TEeff2      = M_prod./(B_prod + MZl_prod + LZl_prod);
    TEeff2(2,:) = L_prod./(B_prod + MZl_prod + LZl_prod);
    
    %%
    TL(e,:) = TLprey(8,:);        %Dem
    TEc(e,:) = TEcon(2,:);        %Large
    TEp1(e,:) = TEprod_grp(3,:);  %Large
    TEp2(e,:) = TEprod(5,:);      %Dem
    TEe1(e,:) = TEeff1(2,:);      %Large
    TEe2(e,:) = TEeff2(2,:);      %Large
    
    save([dpath sname 'locs_' harv '_lastyr_TE_TL.mat'],...
    'Pcon','Fcon','Dcon','Pprod','Fprod','Dprod','B_prod','MZ_prod',...
    'LZ_prod','MZl_prod','LZl_prod','conF','conP','conD','conZm',...
    'conZl','conB','TLprey','Con','TEcon','TEprod_grp','TEprod',...
    'TEeff1','TEeff2');
    
    
end
cfile2 = ['Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
save([datap 'Bio_rates/' cfile2 '_locs_DprefTests_TE_TL.mat']);

%% Plots

for s=1:length(spots)
    f1=figure(1);
    subplot(4,4,s)
    plot(dees,TL(:,s),'b','LineWidth',2)
    %set(gca,'XTick',bees,'XTickLabel',bees,'YTick',cmaxs,'YTickLabel',cmaxs)
    xlabel('D pref')
    ylabel('LD TL')
    title(spots{s})
    
    f2=figure(2);
    subplot(4,4,s)
    plot(dees,TEp2(:,s),'b','LineWidth',2)
    %set(gca,'XTick',bees,'XTickLabel',bees,'YTick',cmaxs,'YTickLabel',cmaxs)
    xlabel('D pref')
    ylabel('LD TE prod')
    title(spots{s})
    
end
print(f1,'-dpng',[figp cfile2 '_Dpref_LD_TL.png'])
print(f2,'-dpng',[figp cfile2 '_Dpref_LD_TEprod.png'])
