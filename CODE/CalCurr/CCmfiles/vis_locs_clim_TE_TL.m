% Visualize output of POEM Climatology at single locations
% 150 years, monthly means saved
% Transfer efficiency ("effective") and trophic level

clear all
close all

warning off

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

dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
sname = 'Clim_';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);
load([dpath sname 'locs_' harv '.mat'])

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

mass = [M_s;M_m;M_l];
mass = repmat(mass,1,length(spots));
L = [L_s;L_m;L_l];

A = 4.39;
fc = 0.2;
f0 = 0.6;
epsassim = 0.7;
n = 3/4;

w = logspace(-3, 5);
AvailEnergy = A*w.^n;
Consumption = A / (epsassim*(f0-fc)) * w.^n;

%Andersen & Beyer mortality rate per year (natural + predation)
%physiol mort * growth constant * M^-0.25
AB = (0.35 .* 4.5 .* mass.^(-0.25)) ./365;

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%% Zoop and det and npp
cpath = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
gpath='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([gpath 'clim_det_biom_Dmeans_Ytot.mat'])
load([gpath 'clim_npp_Dmeans_Ytot.mat'])
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

mnpp = npp_mean_clim(ID)* 1e-3 * 9.0;
tnpp = npp_tot_clim(ID)* 1e-3 * 9.0;

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
mnpp_locs = mnpp(ids);
tnpp_locs = tnpp(ids);

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
Det_prod = mdet_locs';

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
TEprod(3,:) = M_prod2./(Det_prod); %M_prod2./(B_prod);
TEprod(4,:) = L_prod1./M_prod1;
TEprod(5,:) = L_prod2./(M_prod2+Det_prod); %L_prod2./(M_prod2+B_prod);

figure
bar(TEprod)
ylim([0 0.3])

%% Effective TE
TEeff1      = M_prod./(B_prod + MZ_prod + LZ_prod);
TEeff1(2,:) = L_prod./(B_prod + MZ_prod + LZ_prod);

TEeff2      = M_prod./(B_prod + MZl_prod + LZl_prod);
TEeff2(2,:) = L_prod./(B_prod + MZl_prod + LZl_prod);

TEeff3      = M_prod./(Det_prod + MZ_prod + LZ_prod);
TEeff3(2,:) = L_prod./(Det_prod + MZ_prod + LZ_prod);

TEeff4      = M_prod./(Det_prod + MZl_prod + LZl_prod);
TEeff4(2,:) = L_prod./(Det_prod + MZl_prod + LZl_prod);


%%
save([dpath sname 'locs_' harv '_lastyr_TE_TL.mat'],...
    'Pcon','Fcon','Dcon','Pprod','Fprod','Dprod','B_prod','MZ_prod',...
    'LZ_prod','MZl_prod','LZl_prod','conF','conP','conD','conZm',...
    'conZl','conB','TLprey','Con','TEcon','TEprod_grp','TEprod',...
    'TEeff1','TEeff2','TEeff3','TEeff4');

TEM = real(TEeff2(1,:).^(1/3));
TEL = real(TEeff2(2,:).^(1/4));
Tab=table(spots,TEeff2(1,:)',TEeff2(2,:)',TEM',TEL',...
    'VariableNames',{'Site','TEeffM','TEeffL','TEM','TEL'});
writetable(Tab,[dpath 'Locs_TEeffLoss_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'Locs_TEeffLoss_clim_fished_',harv,'_' cfile '.csv'],'Tab');

%writetable(Tab,['Locs_TEeffLoss_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');

%% Figures
%TEcon
figure(1);
subplot(2,1,1)
plot(1:length(spots),TEcon(1,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('mean conTE in final year')
xlabel('Site')
title('Medium/Small')

subplot(2,1,2)
plot(1:length(spots),TEcon(2,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('mean conTE in final year')
xlabel('Site')
title('Large/Medium')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_conTEs.png'])

%% TE prod by habitat and size
figure(2);
subplot(3,2,1)
plot(1:length(spots),TEprod(1,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:2:length(spots),'XTickLabel',spots(1:2:end))
title('S')

subplot(3,2,3)
plot(1:length(spots),TEprod(2,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',2:2:length(spots),'XTickLabel',spots(2:2:end))
ylabel('mean prodTE in final year')
title('MF & MP')

subplot(3,2,4)
plot(1:length(spots),TEprod(3,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',2:2:length(spots),'XTickLabel',spots(2:2:end))
title('MD')

subplot(3,2,5)
plot(1:length(spots),TEprod(4,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:2:length(spots),'XTickLabel',spots(1:2:end))
xlabel('Site')
title('LP')

subplot(3,2,6)
plot(1:length(spots),TEprod(5,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:2:length(spots),'XTickLabel',spots(1:2:end))
xlabel('Site')
title('LD')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_prodTEs.png'])

%% TE eff
figure(3);
subplot(2,1,1)
plot(1:length(spots),TEeff1(1,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Medium')
xlabel('Site')
title('mean TE eff in final year')

subplot(2,1,2)
plot(1:length(spots),TEeff1(2,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Large')
xlabel('Site')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_effTEs.png'])

%Losses
figure(4);
subplot(2,1,1)
plot(1:length(spots),TEeff2(1,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Medium')
xlabel('Site')
title('mean TE eff in final year')

subplot(2,1,2)
plot(1:length(spots),TEeff2(2,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Large')
xlabel('Site')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_effTEs_losses.png'])

%Bio & Det
figure(7);
subplot(2,1,1)
plot(1:length(spots),TEeff3(1,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Medium')
xlabel('Site')
title('mean TE eff in final year')

subplot(2,1,2)
plot(1:length(spots),TEeff3(2,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Large')
xlabel('Site')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_effTEs_det.png'])

%Losses & Det
figure(8);
subplot(2,1,1)
plot(1:length(spots),TEeff4(1,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Medium')
xlabel('Site')
title('mean TE eff in final year')

subplot(2,1,2)
plot(1:length(spots),TEeff4(2,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',spots)
ylabel('Large')
xlabel('Site')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_effTEs_losses_det.png'])

%% TL prey
figure(5);
subplot(3,1,1)
plot(1-0.25:16,TLprey(1,:),'sk',...
    'MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:16,TLprey(2,:),'sk',...
    'MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot(1+0.25:17,TLprey(3,:),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
set(gca,'XTick',1:16,'XTickLabel',spots);
title('S')

subplot(3,1,2)
plot(1-0.25:16,TLprey(4,:),'sk',...
    'MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:16,TLprey(5,:),'sk',...
    'MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot(1+0.25:17,TLprey(6,:),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
set(gca,'XTick',1:16,'XTickLabel',spots);
ylabel('Mean TL in final year')
title('M')

subplot(3,1,3)
plot(1:16,TLprey(7,:),'sk',...
    'MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
plot(1+0.25:17,TLprey(8,:),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
set(gca,'XTick',1:16,'XTickLabel',spots);
title('L')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_TL.png'])

%%
figure(6);
plot(1+0.25:17,TLprey(8,:),'sk',...
    'MarkerFaceColor',cmap_ppt(2,:),...
    'MarkerSize',15); hold on;
xlim([0 17])
set(gca,'XTick',1:16,'XTickLabel',spots);
ylabel('Mean TL in final year')
title('L')
stamp(cfile)
print('-dpng',[fpath sname 'locs_' harv '_LD_TL.png'])




