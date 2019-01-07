% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
%figp = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
sname2 = '';

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

%%
dp1 = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
dp6 = 'Dc_2B_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
dpath1 = [datap char(dp1) '/'];
fpath = figp;
dpath6 = [datap char(dp6) '/'];
cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400_Btest';

load([dpath1 sname 'lastyr_sum_mean_biom.mat'],'all_mean');
fishsp1 = squeeze(nansum(all_mean));
sumspec1 = squeeze(nansum(nansum(all_mean)));
clear all_mean

load([dpath6 sname 'lastyr_sum_mean_biom.mat'],'all_mean');
fishsp6 = squeeze(nansum(all_mean));
sumspec6 = squeeze(nansum(nansum(all_mean)));
clear all_mean

%%
figure(1);
plot((1-0.2):11,log10(fishsp1(1,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:11,log10(fishsp6(1,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
legend('1B','2B')
xlim([0 12])
ylim([-2 2])
set(gca,'XTick',1:11,'XTickLabel',[])
for n=1:11
    text(n,-2.2,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 F Mean Biom (g m^-^2) in final year')
title('All F stages')
stamp(cfile)
print('-dpng',[fpath cfile '_All_oneloc_tot_mean_Fbiomass_1B_2B.png'])

figure(2);
plot((1-0.2):11,log10(fishsp1(2,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:11,log10(fishsp6(2,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
legend('1B','2B')
xlim([0 12])
ylim([-2 2])
set(gca,'XTick',1:11,'XTickLabel',[])
for n=1:11
    text(n,-2.2,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 P Mean Biom (g m^-^2) in final year')
title('All P stages')
stamp(cfile)
print('-dpng',[fpath cfile '_All_oneloc_tot_mean_Dbiomass_1B_2B.png'])

figure(3);
plot((1-0.2):11,log10(fishsp1(3,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:11,log10(fishsp6(3,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
legend('1B','2B')
xlim([0 12])
ylim([-2 2])
set(gca,'XTick',1:11,'XTickLabel',[])
for n=1:11
    text(n,-2.2,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 D Mean Biom (g m^-^2) in final year')
title('All D stages')
stamp(cfile)
print('-dpng',[fpath cfile '_All_oneloc_tot_mean_Dbiomass_1B_2B.png'])

figure(4);
plot((1-0.2):11,log10(fishsp1(4,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1:11,log10(fishsp6(4,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
legend('1B','2B')
xlim([0 12])
ylim([-2 2])
set(gca,'XTick',1:11,'XTickLabel',[])
for n=1:11
    text(n,-2.2,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 B Mean Biom (g m^-^2) in final year')
title('All B sizes')
stamp(cfile)
print('-dpng',[fpath cfile '_All_oneloc_tot_mean_Bbiomass_1B_2B.png'])

%%
figure(5);
plot(1:11,log10(sumspec1),'sk','MarkerFaceColor',cmap_ppt(3,:),...
    'MarkerSize',15); hold on;
plot(1.1:11.1,log10(sumspec6),'sk','MarkerFaceColor',cmap_ppt(1,:),...
    'MarkerSize',15); hold on;
legend('1B','2B')
xlim([0 12])
ylim([-2 2])
set(gca,'XTick',1:11,'XTickLabel',[])
for n=1:11
    text(n,-2.1,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All fishes and stages')
stamp(cfile)
print('-dpng',[fpath cfile '_All_oneloc_tot_mean_biomass_spec_1B_2B.png'])

