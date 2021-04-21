% Calc residence time
% Residence time = 1 / input
% or             = 1 / output
% Total inputs: rec, nu
% Total outputs: gamma, rep, nmort, die (pred), yield (fishing)
% Use per biomass rates (g/g/m2/d)

clear all
close all

%%
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/MIP/NC/Matlab_new_size/param_ensemble/',cfile,'/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught','fmort'};
cols=cols';
spots=spots';
spots{16} = 'PUP';

%%
load([fpath 'Locs_Climatol_All_fish03_means_mzpref_search.mat']);

%% add gains and losses
% +nu is g/g/m2/d; rec is g/m2/d
SF_in = max(0,nSF) + rSF ./ SF;
SP_in = max(0,nSP) + rSP ./ SP;
SD_in = max(0,nSD) + rSD ./ SD;
MF_in = max(0,nMF) + rMF ./ MF;
MP_in = max(0,nMP) + rMP ./ MP;
MD_in = max(0,nMD) + rMD ./ MD;
LP_in = max(0,nLP) + rLP ./ LP;
LD_in = max(0,nLD) + rLD ./ LD;

% (-nu,) gamma, rep, nmort, pred, fmort are g/g/m2/d; mort, die, yieLD are g/m2/d
SF_out = abs(min(0,nSF)) + gSF + mSF + pSF;
SP_out = abs(min(0,nSP)) + gSP + mSP + pSP;
SD_out = abs(min(0,nSD)) + gSD + mSD + pSD;
MF_out = abs(min(0,nMF)) + gMF + mMF + pMF + fMF + eMF;
MP_out = abs(min(0,nMP)) + gMP + mMP + pMP + fMP;
MD_out = abs(min(0,nMD)) + gMD + mMD + pMD + fMD;
LP_out = abs(min(0,nLP)) + gLP + mLP + pLP + fLP + eLP;
LD_out = abs(min(0,nLD)) + gLD + mLD + pLD + fLD + eLD;

%% v1
SF_res1 = 1 ./ SF_in;
SP_res1 = 1 ./ SP_in;
SD_res1 = 1 ./ SD_in;
MF_res1 = 1 ./ MF_in;
MP_res1 = 1 ./ MP_in;
MD_res1 = 1 ./ MD_in;
LP_res1 = 1 ./ LP_in;
LD_res1 = 1 ./ LD_in;

% SF_res1(isinf(SF_res1(:))) = NaN;
% SP_res1(isinf(SP_res1(:))) = NaN;
% SD_res1(isinf(SD_res1(:))) = NaN;
% MF_res1(isinf(MF_res1(:))) = NaN;
% MP_res1(isinf(MP_res1(:))) = NaN;
% MD_res1(isinf(MD_res1(:))) = NaN;
% LP_res1(isinf(LP_res1(:))) = NaN;
% LD_res1(isinf(LD_res1(:))) = NaN;

%% v2
SF_res2 = 1 ./ SF_out;
SP_res2 = 1 ./ SP_out;
SD_res2 = 1 ./ SD_out;
MF_res2 = 1 ./ MF_out;
MP_res2 = 1 ./ MP_out;
MD_res2 = 1 ./ MD_out;
LP_res2 = 1 ./ LP_out;
LD_res2 = 1 ./ LD_out;

%% Save
save([fpath 'Residence_time_means_Climatol_' harv '_locs_mzpref.mat'],...
    'SF_res1','SP_res1','SD_res1','MF_res1','MP_res1','MD_res1','LP_res1','LD_res1',...
    'SF_res2','SP_res2','SD_res2','MF_res2','MP_res2','MD_res2','LP_res2','LD_res2')

%% Histograms
edges = [0:30:360 547 730 912 1095];
edge2 = -4:2;

%biomass
figure(10)
subplot(3,3,1)
histogram(log10(SF),edge2)
title('SF')

subplot(3,3,2)
histogram(log10(SP),edge2)
title('SP')

subplot(3,3,3)
histogram(log10(SD),edge2)
title('SD')

subplot(3,3,4)
histogram(log10(MF),edge2)
title('MF')

subplot(3,3,5)
histogram(log10(MP),edge2)
title('MP')

subplot(3,3,6)
histogram(log10(MD),edge2)
title('MD')

subplot(3,3,8)
histogram(log10(LP),edge2)
title('LP')

subplot(3,3,9)
histogram(log10(LD),edge2)
title('LD')

%res 1
figure(11)
subplot(3,3,1)
histogram((SF_res1),edges)
title('SF')

subplot(3,3,2)
histogram((SP_res1),edges)
title('SP')

subplot(3,3,3)
histogram((SD_res1),edges)
title('SD')

subplot(3,3,4)
histogram((MF_res1),edges)
title('MF')

subplot(3,3,5)
histogram((MP_res1),edges)
title('MP')

subplot(3,3,6)
histogram((MD_res1),edges)
title('MD')

subplot(3,3,8)
histogram((LP_res1),edges)
title('LP')

subplot(3,3,9)
histogram((LD_res1),edges)
title('LD')

%res 2
figure(12)
subplot(3,3,1)
histogram((SF_res2),edges)
title('SF')

subplot(3,3,2)
histogram((SP_res2),edges)
title('SP')

subplot(3,3,3)
histogram((SD_res2),edges)
title('SD')

subplot(3,3,4)
histogram((MF_res2),edges)
title('MF')

subplot(3,3,5)
histogram((MP_res2),edges)
title('MP')

subplot(3,3,6)
histogram((MD_res2),edges)
title('MD')

subplot(3,3,8)
histogram((LP_res2),edges)
title('LP')

subplot(3,3,9)
histogram((LD_res2),edges)
title('LD')

%% Figures
nk=length(mzp);

Fbio = SF + MF;
Pbio = SP + MP + LP;
Dbio = SD + MD + LD;

cm=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255];
set(groot,'defaultAxesColorOrder',cm);

%0.75 0.75 0.75;... %lt grey
%1 1 0;...      %yellow

%% Only use 3 domain examples: EBS(10), PUp(16), HOT(13)
% sid = [10;16;13];
% for s=1:length(sid)
%     domain = sid(s);
%     loc = spots{domain};
%     lname = [loc '_'];

% Biomass
f1=figure(1);
plot(mzp,log10(Fbio),'LineWidth',2)
%yaxis([-4 2])
legend(spots)
legend('location','eastoutside')
title('F biomass')

f2=figure(2);
plot(mzp,log10(Pbio),'LineWidth',2)
%yaxis([-4 2])
legend(spots)
legend('location','eastoutside')
title('P biomass')

f3=figure(3);
plot(mzp,log10(Dbio),'LineWidth',2)
%yaxis([-4 2])
legend(spots)
legend('location','eastoutside')
title('D biomass')

%% Small res
f4=figure(4);
plot(mzp,SF_res2,'LineWidth',2)
ylim([1 100])
legend(spots)
legend('location','eastoutside')
title('SF res')

% f5=figure(5);
% plot(mzp,SP_res2,'LineWidth',2)
% ylim([1 100])
% legend(spots)
% legend('location','eastoutside')
% title('SP res')
% 
% f6=figure(6);
% plot(mzp,SD_res2,'LineWidth',2)
% ylim([1 100])
% legend(spots)
% legend('location','eastoutside')
% title('SD res')

%% Medium res
f7=figure(7);
subplot(2,2,1)
plot(mzp,MF_res2,'LineWidth',2)
ylim([10 300])
xlim([0.1 0.95])
title('MF res')

subplot(2,2,2)
plot(mzp,MP_res2,'LineWidth',2)
ylim([10 300])
xlim([0.1 0.95])
title('MP res')

subplot(2,1,2)
plot(mzp,MD_res2,'LineWidth',2)
ylim([10 600])
xlim([0.1 0.95])
legend(spots)
legend('location','eastoutside')
title('MD res')


%% Large res
f8=figure(8);
subplot(2,1,1)
plot(mzp,LP_res2,'LineWidth',2)
%ylim([10 600])
xlim([0.1 0.95])
legend(spots)
legend('location','eastoutside')
title('LP res')

subplot(2,1,2)
plot(mzp,LD_res2,'LineWidth',2)
%ylim([10 600])
xlim([0.1 0.95])
legend(spots)
legend('location','eastoutside')
title('LD res')

%%
print(f1,'-dpng',[ppath 'Fbiomass_locs_mzpref.png'])
print(f2,'-dpng',[ppath 'Pbiomass_locs_mzpref.png'])
print(f3,'-dpng',[ppath 'Dbiomass_locs_mzpref.png'])
print(f4,'-dpng',[ppath 'SFres2_locs_mzpref.png'])
% print(f5,'-dpng',[ppath 'SPres2_locs_mzpref.png'])
% print(f6,'-dpng',[ppath 'SDres2_locs_mzpref.png'])
print(f7,'-dpng',[ppath 'Mres2_locs_mzpref.png'])
print(f8,'-dpng',[ppath 'Lres2_locs_mzpref.png'])



