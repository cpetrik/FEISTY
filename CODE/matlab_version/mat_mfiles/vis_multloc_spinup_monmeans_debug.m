% Visualize output of POEM
% Spinup at 100 locations
% 50 years
% Saved as mat files

% clear all
% close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Mat_runs/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC275_RE0500';

fpath=['/Volumes/GFDL/NC/Matlab_runs/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_spinup_' cfile '.mat']);

%% colors
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
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

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in time
y = time;

% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup Pelagic Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_P_time.png'])

%% Planktivore
sf_tmean=sf_tmean(1:length(y));
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup Forage Fishes')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
legend('Immature','Adults')
legend('location','southeast')
stamp(cfile)

subplot(3,1,2)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Immature')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,1,3)
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_F_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup Demersal Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_D_time.png'])

%% All size classes of all

figure(5)
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup')
stamp(cfile)
print('-dpng',[ppath 'Spinup_all_sizes.png'])


%% Plots in space
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

%% DEBUGGING -------------------------------------------------------------
% Bio
for n=10:12
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.bio(:,n),X,Y);
    Zsf=griddata(grid(:,2),grid(:,3),SF.bio(:,n),X,Y);
    Zsd=griddata(grid(:,2),grid(:,3),SD.bio(:,n),X,Y);
    Zmp=griddata(grid(:,2),grid(:,3),MP.bio(:,n),X,Y);
    Zmf=griddata(grid(:,2),grid(:,3),MF.bio(:,n),X,Y);
    Zmd=griddata(grid(:,2),grid(:,3),MD.bio(:,n),X,Y);
    Zlp=griddata(grid(:,2),grid(:,3),LP.bio(:,n),X,Y);
    Zld=griddata(grid(:,2),grid(:,3),LD.bio(:,n),X,Y);
    Zb=griddata(grid(:,2),grid(:,3),BENT.bio(:,n),X,Y);
    
    figure(50)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zb)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean benthic biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_BENT' num2str(n) '.png'])
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsp)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SP' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsf)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SF' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsd)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SD' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmp)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MP' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmf)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
    title(['log10 mean Adult F biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MF' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmd)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MD' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zlp)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
    title(['log10 mean Adult P biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LP' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zld)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
    title(['log10 mean Adult D biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LD' num2str(n) '.png'])
end

%% Con
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

mass = [M_s;M_m;M_l];

A = 4.39;
fc = 0.2;
f0 = 0.6;
epsassim = 0.6;
n = 3/4;

w = logspace(-3, 5);

AvailEnergy = A*w.^n;
Consumption = A / (epsassim*(f0-fc)) * w.^n;
Consumption2 = A / (epsassim*(f0-fc)) * mass.^n;
lgCon2 = log10(Consumption2);

% Consump g/g/d --> g/d --> g/y

for n=10:12
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.con(:,n),X,Y) .* M_s .* 365;
    Zsf=griddata(grid(:,2),grid(:,3),SF.con(:,n),X,Y) .* M_s .* 365;
    Zsd=griddata(grid(:,2),grid(:,3),SD.con(:,n),X,Y) .* M_s .* 365;
    Zmp=griddata(grid(:,2),grid(:,3),MP.con(:,n),X,Y) .* M_m .* 365;
    Zmf=griddata(grid(:,2),grid(:,3),MF.con(:,n),X,Y) .* M_m .* 365;
    Zmd=griddata(grid(:,2),grid(:,3),MD.con(:,n),X,Y) .* M_m .* 365;
    Zlp=griddata(grid(:,2),grid(:,3),LP.con(:,n),X,Y) .* M_l .* 365;
    Zld=griddata(grid(:,2),grid(:,3),LD.con(:,n),X,Y) .* M_l .* 365;
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsp)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SPcon' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsf)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SFcon' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsd)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SDcon' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmp)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPcon' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmf)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    title(['log10 mean Adult F con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MFcon' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmd)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDcon' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zlp)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult P con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LPcon' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zld)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult D con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LDcon' num2str(n) '.png'])
end

%% Clev
%Time
clev_tmean(1,:)=mean(SF.clev(1:600),1);
clev_tmean(2,:)=mean(SP.clev,1);
clev_tmean(3,:)=mean(SD.clev,1);
clev_tmean(4,:)=mean(MF.clev,1);
clev_tmean(5,:)=mean(MP.clev,1);
clev_tmean(6,:)=mean(MD.clev,1);
clev_tmean(7,:)=mean(LP.clev,1);
clev_tmean(8,:)=mean(LD.clev,1);

figure
plot(1:600,clev_tmean)
legend('SF','SP','SD','MF','MP','MD','LP','LD')
legend('location','eastoutside')
xlim([1 20])
xlabel('Time (mo)')
ylim([0 1.1])
ylabel('C/Cmax')
title('Spinup')
stamp(cfile)
print('-dpng',[ppath 'Spinup_all_sizes_clev.png'])

% Last year
lyr=time((end-12+1):end);
sp_cmean=mean(SP.clev(:,lyr),2);
sf_cmean=mean(SF.clev(:,lyr),2);
sd_cmean=mean(SD.clev(:,lyr),2);
mp_cmean=mean(MP.clev(:,lyr),2);
mf_cmean=mean(MF.clev(:,lyr),2);
md_cmean=mean(MD.clev(:,lyr),2);
lp_cmean=mean(LP.clev(:,lyr),2);
ld_cmean=mean(LD.clev(:,lyr),2);

Zsp=griddata(grid(:,2),grid(:,3),sp_cmean,X,Y) ;
Zsf=griddata(grid(:,2),grid(:,3),sf_cmean,X,Y);
Zsd=griddata(grid(:,2),grid(:,3),sd_cmean,X,Y);
Zmp=griddata(grid(:,2),grid(:,3),mp_cmean,X,Y);
Zmf=griddata(grid(:,2),grid(:,3),mf_cmean,X,Y);
Zmd=griddata(grid(:,2),grid(:,3),md_cmean,X,Y);
Zlp=griddata(grid(:,2),grid(:,3),lp_cmean,X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_cmean,X,Y);

% sp
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zsp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean Larval P C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SPclev.png'])

% sf
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zsf); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval F C/Cmax')
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SFclev.png'])

% sd
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zsd); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean Larval D C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SDclev.png'])

% mp
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zmp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean Juvenile P C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MPclev.png'])

% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zmf); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_plot(-111,5,'o','w','MarkerSize',10);
title(['mean Adult F C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MFclev.png'])

% md
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zmd); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['Juvenile D C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MDclev.png'])

% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zlp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean Adult P C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LPclev.png'])

% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zld); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean Adult D C/Cmax'])
colormap('jet')
colorbar('h')
%caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LDclev.png'])

%% Zoop
load('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2005.mat');

% Last year
zm_mean=mean(COBALT.Zm,2);
zl_mean=mean(COBALT.Zl,2);
dzm_mean=mean(COBALT.dZm,2);
dzl_mean=mean(COBALT.dZl,2);

Zm=griddata(grid(:,2),grid(:,3),zm_mean,X,Y) ;
Zl=griddata(grid(:,2),grid(:,3),zl_mean,X,Y);
dZm=griddata(grid(:,2),grid(:,3),dzm_mean,X,Y) ;
dZl=griddata(grid(:,2),grid(:,3),dzl_mean,X,Y);

% mz
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zm); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean MZ'])
colormap('jet')
colorbar('h')
%caxis([0 10])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MZ.png'])

% lz
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zl); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean LZ')
colormap('jet')
colorbar('h')
%caxis([0 10])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LZ.png'])

% mz
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,dZm); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title(['mean dMZ'])
colormap('jet')
colorbar('h')
%caxis([0 10])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_dMZ.png'])

% lz
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,dZl); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean dLZ')
colormap('jet')
colorbar('h')
%caxis([0 10])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_dLZ.png'])



%% Gamma/rep
for n=10:12
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.gamma(:,n),X,Y);
    Zsf=griddata(grid(:,2),grid(:,3),SF.gamma(:,n),X,Y);
    Zsd=griddata(grid(:,2),grid(:,3),SD.gamma(:,n),X,Y);
    Zmp=griddata(grid(:,2),grid(:,3),MP.gamma(:,n),X,Y);
    Zmf=griddata(grid(:,2),grid(:,3),MF.rep(:,n),X,Y);
    Zmd=griddata(grid(:,2),grid(:,3),MD.gamma(:,n),X,Y);
    Zlp=griddata(grid(:,2),grid(:,3),LP.rep(:,n),X,Y);
    Zld=griddata(grid(:,2),grid(:,3),LD.rep(:,n),X,Y);
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsp+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SPgamma' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsf+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SFgamma' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zsd+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SDgamma' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmp+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPgamma' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmf+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    title(['log10 mean Adult F rep (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MFrep' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zmd+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDgamma' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zlp+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult P rep (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LPrep' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log10(Zld+eps)); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult D rep (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LDrep' num2str(n) '.png'])
end

%% Gamma/rep
for n=10:12
    close all
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['mean Juvenile P gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-1e-25 0])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPgamma' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmd); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['mean Juvenile D gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-1e-25 0])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDgamma' num2str(n) '.png'])
    
end
%%
figure
plot(SP.gamma(:,12),'.')
test=log10(SP.gamma(:,12));
figure
plot(test,'.')

figure
plot(MP.gamma(:,12),'.')
test=log10(MP.gamma(:,12));
figure
plot(test,'.')

%% Nu
for n=10:12
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.nu(:,n),X,Y);
    Zsf=griddata(grid(:,2),grid(:,3),SF.nu(:,n),X,Y);
    Zsd=griddata(grid(:,2),grid(:,3),SD.nu(:,n),X,Y);
    Zmp=griddata(grid(:,2),grid(:,3),MP.nu(:,n),X,Y);
    Zmf=griddata(grid(:,2),grid(:,3),MF.nu(:,n),X,Y);
    Zmd=griddata(grid(:,2),grid(:,3),MD.nu(:,n),X,Y);
    Zlp=griddata(grid(:,2),grid(:,3),LP.nu(:,n),X,Y);
    Zld=griddata(grid(:,2),grid(:,3),LD.nu(:,n),X,Y);
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zsp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SPnu' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zsf); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SFnu' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zsd); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SDnu' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPnu' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmf); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    title(['log10 mean Adult F nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MFnu' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmd); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDnu' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zlp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult P nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LPnu' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zld); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult D nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LDnu' num2str(n) '.png'])
end

%% Find negative values
mos=10:12;
tend=length(mos);
bio=NaN*ones(9,tend);
con=NaN*ones(8,tend);
die=NaN*ones(8,tend);
gamma=NaN*ones(8,tend);
rec=NaN*ones(8,tend);

for k=1:tend
    n=mos(k);
    %bio
    id = (SP.bio(:,n) < 0);
    bio(1,k) = sum(id);
    id = (SF.bio(:,n) < 0);
    bio(2,k) = sum(id);
    id = (SD.bio(:,n) < 0);
    bio(3,k) = sum(id);
    id = (MP.bio(:,n) < 0);
    bio(4,k) = sum(id);
    id = (MF.bio(:,n) < 0);
    bio(5,k) = sum(id);
    id = (MD.bio(:,n) < 0);
    bio(6,k) = sum(id);
    id = (LP.bio(:,n) < 0);
    bio(7,k) = sum(id);
    id = (LD.bio(:,n) < 0);
    bio(8,k) = sum(id);
    id = (BENT.bio(:,n) < 0);
    bio(9,k) = sum(id);
    %con
    id = (SP.con(:,n) < 0);
    con(1,k) = sum(id);
    id = (SF.con(:,n) < 0);
    con(2,k) = sum(id);
    id = (SD.con(:,n) < 0);
    con(3,k) = sum(id);
    id = (MP.con(:,n) < 0);
    con(4,k) = sum(id);
    id = (MF.con(:,n) < 0);
    con(5,k) = sum(id);
    id = (MD.con(:,n) < 0);
    con(6,k) = sum(id);
    id = (LP.con(:,n) < 0);
    con(7,k) = sum(id);
    id = (LD.con(:,n) < 0);
    con(8,k) = sum(id);
    %die
    id = (SP.die(:,n) < 0);
    die(1,k) = sum(id);
    id = (SF.die(:,n) < 0);
    die(2,k) = sum(id);
    id = (SD.die(:,n) < 0);
    die(3,k) = sum(id);
    id = (MP.die(:,n) < 0);
    die(4,k) = sum(id);
    id = (MF.die(:,n) < 0);
    die(5,k) = sum(id);
    id = (MD.die(:,n) < 0);
    die(6,k) = sum(id);
    id = (LP.die(:,n) < 0);
    die(7,k) = sum(id);
    id = (LD.die(:,n) < 0);
    die(8,k) = sum(id);
    %gamma/rep
    id = (SP.gamma(:,n) < 0);
    gamma(1,k) = sum(id);
    id = (SF.gamma(:,n) < 0);
    gamma(2,k) = sum(id);
    id = (SD.gamma(:,n) < 0);
    gamma(3,k) = sum(id);
    id = (MP.gamma(:,n) < 0);
    gamma(4,k) = sum(id);
    id = (MF.rep(:,n) < 0);
    gamma(5,k) = sum(id);
    id = (MD.gamma(:,n) < 0);
    gamma(6,k) = sum(id);
    id = (LP.rep(:,n) < 0);
    gamma(7,k) = sum(id);
    id = (LD.rep(:,n) < 0);
    gamma(8,k) = sum(id);
    %rec
    id = (SP.rec(:,n) < 0);
    rec(1,k) = sum(id);
    id = (SF.rec(:,n) < 0);
    rec(2,k) = sum(id);
    id = (SD.rec(:,n) < 0);
    rec(3,k) = sum(id);
    id = (MP.rec(:,n) < 0);
    rec(4,k) = sum(id);
    id = (MF.rec(:,n) < 0);
    rec(5,k) = sum(id);
    id = (MD.rec(:,n) < 0);
    rec(6,k) = sum(id);
    id = (LP.rec(:,n) < 0);
    rec(7,k) = sum(id);
    id = (LD.rec(:,n) < 0);
    rec(8,k) = sum(id);
    
end
%%
figure
plot(mos,bio)
%ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD','B')
legend('location','northwest')
print('-dpng',[ppath 'neg_biotest.png'])

figure
plot(mos,con)
%ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'neg_contest.png'])

figure
plot(mos,die)
%ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'neg_dietest.png'])

figure
plot(mos,gamma)
%ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'neg_gammatest.png'])

figure
plot(mos,rec)
%ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'neg_rectest.png'])

%% Find NaNs
mos=10:12;
tend=length(mos);
bio=NaN*ones(9,tend);
con=NaN*ones(8,tend);
die=NaN*ones(8,tend);
gamma=NaN*ones(8,tend);
rec=NaN*ones(8,tend);
nu=NaN*ones(8,tend);

for k=1:tend
    n=mos(k);
    %bio
    id = isnan(SP.bio(:,n));
    bio(1,k) = sum(id);
    id = isnan(SF.bio(:,n));
    bio(2,k) = sum(id);
    id = isnan(SD.bio(:,n));
    bio(3,k) = sum(id);
    id = isnan(MP.bio(:,n));
    bio(4,k) = sum(id);
    id = isnan(MF.bio(:,n));
    bio(5,k) = sum(id);
    id = isnan(MD.bio(:,n));
    bio(6,k) = sum(id);
    id = isnan(LP.bio(:,n));
    bio(7,k) = sum(id);
    id = isnan(LD.bio(:,n));
    bio(8,k) = sum(id);
    id = isnan(BENT.bio(:,n));
    bio(9,k) = sum(id);
    %con
    id = isnan(SP.con(:,n));
    con(1,k) = sum(id);
    id = isnan(SF.con(:,n));
    con(2,k) = sum(id);
    id = isnan(SD.con(:,n));
    con(3,k) = sum(id);
    id = isnan(MP.con(:,n));
    con(4,k) = sum(id);
    id = isnan(MF.con(:,n));
    con(5,k) = sum(id);
    id = isnan(MD.con(:,n));
    con(6,k) = sum(id);
    id = isnan(LP.con(:,n));
    con(7,k) = sum(id);
    id = isnan(LD.con(:,n));
    con(8,k) = sum(id);
    %die
    id = isnan(SP.die(:,n));
    die(1,k) = sum(id);
    id = isnan(SF.die(:,n));
    die(2,k) = sum(id);
    id = isnan(SD.die(:,n));
    die(3,k) = sum(id);
    id = isnan(MP.die(:,n));
    die(4,k) = sum(id);
    id = isnan(MF.die(:,n));
    die(5,k) = sum(id);
    id = isnan(MD.die(:,n));
    die(6,k) = sum(id);
    id = isnan(LP.die(:,n));
    die(7,k) = sum(id);
    id = isnan(LD.die(:,n));
    die(8,k) = sum(id);
    %gamma/rep
    id = isnan(SP.gamma(:,n));
    gamma(1,k) = sum(id);
    id = isnan(SF.gamma(:,n));
    gamma(2,k) = sum(id);
    id = isnan(SD.gamma(:,n));
    gamma(3,k) = sum(id);
    id = isnan(MP.gamma(:,n));
    gamma(4,k) = sum(id);
    id = isnan(MF.rep(:,n));
    gamma(5,k) = sum(id);
    id = isnan(MD.gamma(:,n));
    gamma(6,k) = sum(id);
    id = isnan(LP.rep(:,n));
    gamma(7,k) = sum(id);
    id = isnan(LD.rep(:,n));
    gamma(8,k) = sum(id);
    %rec
    id = isnan(SP.rec(:,n));
    rec(1,k) = sum(id);
    id = isnan(SF.rec(:,n));
    rec(2,k) = sum(id);
    id = isnan(SD.rec(:,n));
    rec(3,k) = sum(id);
    id = isnan(MP.rec(:,n));
    rec(4,k) = sum(id);
    id = isnan(MF.rec(:,n));
    rec(5,k) = sum(id);
    id = isnan(MD.rec(:,n));
    rec(6,k) = sum(id);
    id = isnan(LP.rec(:,n));
    rec(7,k) = sum(id);
    id = isnan(LD.rec(:,n));
    rec(8,k) = sum(id);
    %nu
    id = isnan(SP.nu(:,n));
    nu(1,k) = sum(id);
    id = isnan(SF.nu(:,n));
    nu(2,k) = sum(id);
    id = isnan(SD.nu(:,n));
    nu(3,k) = sum(id);
    id = isnan(MP.nu(:,n));
    nu(4,k) = sum(id);
    id = isnan(MF.nu(:,n));
    nu(5,k) = sum(id);
    id = isnan(MD.nu(:,n));
    nu(6,k) = sum(id);
    id = isnan(LP.nu(:,n));
    nu(7,k) = sum(id);
    id = isnan(LD.nu(:,n));
    nu(8,k) = sum(id);
    
end
%%
figure
bar(bio')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD','B')
legend('location','northwest')
print('-dpng',[ppath 'nan_biotest.png'])

figure
bar(con')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[ppath 'nan_contest.png'])

figure
bar(die')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[ppath 'nan_dietest.png'])

figure
bar(gamma')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[ppath 'nan_gammatest.png'])

figure
bar(rec')
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[ppath 'nan_rectest.png'])

figure
bar(nu')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[ppath 'nan_nutest.png'])

%%
% n=253;
% id=find(isnan(SP.bio(:,n)));    %43320
% grid(id,2:3);                   %-15.5000lon   60.5000lat
% id2=find(isnan(SP.die(:,n)));
% grid(id2,2:3);
% id3=find(isnan(SP.con(:,n)));
% grid(id3,2:3);
% id4=find(isnan(MP.con(:,n)));
% grid(id4,2:3);
% id5=find(isnan(SP.rec(:,n)))
% grid(id5,2:3)

n=12;
id=find(isnan(SP.bio(:,n)));    
grid(id,2:3);                   
id2=find(isnan(SP.die(:,n)));
grid(id2,2:3);
id3=find(isnan(SP.con(:,n)));   %empty
grid(id3,2:3);
id4=find(isnan(MP.con(:,n)));
grid(id4,2:3);
id5=find(isnan(SP.rec(:,n)))
grid(id5,2:3)

% id5 = 
%        22448
%        33905
%        33906
%        37363
% ans = 
%   129.5000   -6.4424 N Aus
%   -63.5000   14.4587 Caribbean
%   -62.5000   14.4587 Caribbean
%    55.5000   26.5198 Persian Gulf


% END DEBUGGING -------------------------------------------------------------

%%
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');

[ni,nj]=size(geolon_t);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

ocean=NaN*ones(ni,nj);
ocean(grid(:,1))=ones(size(sf_mean));

Zsf(grid(:,1))=sf_mean;
Zsp(grid(:,1))=sp_mean;
Zsd(grid(:,1))=sd_mean;
Zmf(grid(:,1))=mf_mean;
Zmp(grid(:,1))=mp_mean;
Zmd(grid(:,1))=md_mean;
Zlp(grid(:,1))=lp_mean;
Zld(grid(:,1))=ld_mean;
Zb(grid(:,1))=b_mean;

%% ocean cells
figure(55)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,ocean); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('Water cells')
colormap('jet')
colorbar('h')
caxis([1 2])
stamp(cfile)
print('-dpng',[ppath 'Ocean_cells.png'])

% bent
figure(50)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT.png'])

%
mgZb = (Zb/9)*1e3;
figure(51)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(mgZb)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean benthic biomass (mg C m^-^2)')
colormap('jet')
colorbar('h')
caxis([-0.8 2.3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT_mgC.png'])

% sp
figure(11)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SP.png'])

% sf
figure(12)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SF.png'])

% sd
figure(13)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SD.png'])

% mp
figure(14)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MP.png'])

% mf
figure(15)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MF.png'])

% md
figure(16)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MD.png'])

% lp
figure(17)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LP.png'])

% ld
figure(18)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LD.png'])

%% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
FracPDs = Zsp ./ (Zsp+Zsd);
FracPDm = Zmp ./ (Zmp+Zmd);
FracPDl = Zlp ./ (Zlp+Zld);
FracPFs = Zsp ./ (Zsp+Zsf);
FracPFm = Zmp ./ (Zmp+Zmf);
FracPFvDs = (Zsp+Zsf) ./ (Zsp+Zsf+Zsd);
FracPFvDm = (Zmp+Zmf) ./ (Zmp+Zmf+Zmd);

%% ALL
figure(21)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(All)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_All.png'])

% all F
figure(22)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(AllF)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllF.png'])

% all D
figure(23)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(AllD)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllD.png'])

% All P
figure(24)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(AllP)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllP.png'])

% FracPD
figure(25)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPD); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('P:D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPD.png'])

% FracPF
figure(26)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPF); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('P:F mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPF.png'])

%% FracPFvD
figure(27)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFvD); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('(P+F):D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvD.png'])

% FracPDs
figure(28)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPDs); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('SP:SD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPDs.png'])

% FracPFs
figure(29)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFs); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('SP:SF mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFs.png'])

% FracPFvDs
figure(30)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFvDs); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('(SP+SF):SD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvDs.png'])

% FracPDm
figure(31)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPDm); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('MP:MD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPDm.png'])

% FracPFm
figure(32)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFm); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('MP:MF mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFm.png'])

% FracPFvDm
figure(33)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFvDm); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('(MP+MF):MD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvDm.png'])

