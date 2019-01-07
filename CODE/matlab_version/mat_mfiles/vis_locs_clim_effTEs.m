% Visualize output of POEM Climatology at single locations
% 150 years, monthly means saved
% Transfer efficiency ("effective") and trophic level

clear all
close all

warning off

datap = '/Volumes/GFDL/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
BE = 0.075;
sname = 'Climatol_';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);
load([dpath sname harv '_locs.mat'])

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
NPP = mnpp_locs';

S_prod=(SP_prod+SF_prod+SD_prod);
M_prod=(MP_prod+MF_prod+MD_prod);
L_prod=(LP_prod+LD_prod);

Prod = [MZl_prod;LZl_prod;B_prod;S_prod;M_prod;L_prod];



%% Effective TE
TEeff      = M_prod./(B_prod + MZl_prod + LZl_prod);
TEeff(2,:) = L_prod./(B_prod + MZl_prod + LZl_prod);
TEeff(3,:) = L_prod./NPP;
TEeff(4,:) = (B_prod + MZl_prod + LZl_prod)./NPP;

TEeff(5,:) = M_prod./(BE*Det_prod + MZl_prod + LZl_prod);
TEeff(6,:) = L_prod./(BE*Det_prod + MZl_prod + LZl_prod);
TEeff(7,:) = L_prod./NPP;
TEeff(8,:) = (BE*Det_prod + MZl_prod + LZl_prod)./NPP;

%%
save([dpath sname harv '_locs_lastyr_effTEs.mat'],...
    'Pprod','Fprod','Dprod','B_prod','MZ_prod',...
    'LZ_prod','MZl_prod','LZl_prod','TEeff');

TE(1,:) = real(TEeff(1,:).^(1/2));
TE(2,:) = real(TEeff(2,:).^(1/3));
TE(3,:) = real(TEeff(3,:).^(1/4));
TE(4,:) = real(TEeff(4,:));
TE(5,:) = real(TEeff(5,:).^(1/2));
TE(6,:) = real(TEeff(6,:).^(1/3));
TE(7,:) = real(TEeff(7,:).^(1/4));
TE(8,:) = real(TEeff(8,:));
Tab=table(spots,TEeff(1,:)',TEeff(2,:)',TEeff(3,:)',TEeff(4,:)',...
    TEeff(5,:)',TEeff(6,:)',TEeff(7,:)',TEeff(8,:)',...
    'VariableNames',{'Site','TEeff_Mb','TEeff_HTLb','TEeff_Lb','TEeff_LTLb',...
    'TEeff_Md','TEeff_HTLd','TEeff_Ld','TEeff_LTLd'});
writetable(Tab,[dpath 'Locs_TEeff_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'Locs_TEeff_clim_fished_',harv,'_' cfile '.mat'],'Tab');

Tab2=table(spots,TE(1,:)',TE(2,:)',TE(3,:)',TE(4,:)',...
    TE(5,:)',TE(6,:)',TE(7,:)',TE(8,:)',...
    'VariableNames',{'Site','TE_Mb','TE_HTLb','TE_Lb','TE_LTLb',...
    'TE_Md','TE_HTLd','TE_Ld','TE_LTLd'});
writetable(Tab2,[dpath 'Locs_TE_clim_fished_',harv,'_' cfile '.csv'],'Delimiter',',');
save([dpath 'Locs_TE_clim_fished_',harv,'_' cfile '.mat'],'Tab2');


%% Figures
% TE eff
figure(1);
subplot(2,2,1)
plot(1:length(spots),TEeff(1,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TEeff(5,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.01,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('Medium')
ylabel('mean TE eff')

subplot(2,2,2)
plot(1:length(spots),TEeff(2,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TEeff(6,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.001,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('HTL')

subplot(2,2,3)
plot(1:length(spots),TEeff(4,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TEeff(8,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.01,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('LTL')
ylabel('mean TE eff')

subplot(2,2,4)
plot(1:length(spots),TEeff(3,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TEeff(7,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.0001,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('L')
legend('bent','det')
stamp(cfile)
print('-dpng',[fpath sname harv '_locs_effTEs.png'])

%% TE eff converted
figure(2);
subplot(2,2,1)
plot(1:length(spots),TE(1,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TE(5,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.01,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('Medium')
ylabel('mean TE ')

subplot(2,2,2)
plot(1:length(spots),TE(2,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TE(6,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.01,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('HTL')

subplot(2,2,3)
plot(1:length(spots),TE(4,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TE(8,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,-0.01,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('LTL')
ylabel('mean TE ')
legend('bent','det')

subplot(2,2,4)
plot(1:length(spots),TE(3,:),'.b','MarkerSize',25); hold on;
plot(1:length(spots),TE(7,:),'.k','MarkerSize',25); hold on;
xlim([0 length(spots)+1])
set(gca,'XTick',1:length(spots),'XTickLabel',[])
for n=1:16
    text(n,0.04,spots{n},'HorizontalAlignment','right','Rotation',90)
end
title('L')
stamp(cfile)
print('-dpng',[fpath sname harv '_locs_TEs.png'])

