% Calc LME biomass of FEISTY
% CORE-forced
% Observed effort
% Saved as mat files

clear all
close all

%% CORE-forced
load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

ID = GRD.ID;

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90;
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% ESM2M = same grid as CORE
gpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/FEISTY_other/cobalt_data/';
load([gpath 'hindcast_gridspec.mat'],'AREA_OCN');
load([gpath 'lme_mask_esm2m.mat']);
load([cpath 'LME_hist9095_temp_zoop_det.mat'],'lme_ptemp','lme_area');

AREA_OCN = AREA_OCN*510072000*1e6;
AREA_OCN = max(AREA_OCN,1);
tlme = lme_mask_esm2m';

%% fish
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/CORE/'];
ppath = [pp cfile '/CORE/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

harv = 'All_fish_obs'; %All_fish_obs
tharv = 'Observed effort';

load([fpath 'Means_core_',harv,'_' cfile '.mat']);

[nid,nt] = size(ld_mean);

%% Loop over every year to have ts for comp with fishing
lme_mcatch = NaN*ones(66,nt,5);
lme_mbio = NaN*ones(66,nt,9);
lme_sbio = NaN*ones(66,nt,9);

for t=1:nt
    Zsf=NaN*ones(ni,nj);
    Zsp=NaN*ones(ni,nj);
    Zsd=NaN*ones(ni,nj);
    Zmf=NaN*ones(ni,nj);
    Zmp=NaN*ones(ni,nj);
    Zmd=NaN*ones(ni,nj);
    Zlp=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);
    
    Cmf=NaN*ones(ni,nj);
    Cmp=NaN*ones(ni,nj);
    Cmd=NaN*ones(ni,nj);
    Clp=NaN*ones(ni,nj);
    Cld=NaN*ones(ni,nj);
    
    Zsf(ID)=sf_mean(:,t);
    Zsp(ID)=sp_mean(:,t);
    Zsd(ID)=sd_mean(:,t);
    Zmf(ID)=mf_mean(:,t);
    Zmp(ID)=mp_mean(:,t);
    Zmd(ID)=md_mean(:,t);
    Zlp(ID)=lp_mean(:,t);
    Zld(ID)=ld_mean(:,t);
    Zb(ID)=b_mean(:,t);
    
    Cmf(ID)=mf_my(:,t);
    Cmp(ID)=mp_my(:,t);
    Cmd(ID)=md_my(:,t);
    Clp(ID)=lp_my(:,t);
    Cld(ID)=ld_my(:,t);
    
    % g/m2/d --> total g
    Amf_mcatch = Cmf .* AREA_OCN * 365; %mean fish catch per yr
    Amp_mcatch = Cmp .* AREA_OCN * 365;
    Amd_mcatch = Cmd .* AREA_OCN * 365;
    Alp_mcatch = Clp .* AREA_OCN * 365;
    Ald_mcatch = Cld .* AREA_OCN * 365;
    % g/m2 --> total g
    Asf_mean = Zsf .* AREA_OCN;
    Asp_mean = Zsp .* AREA_OCN;
    Asd_mean = Zsd .* AREA_OCN;
    Amf_mean = Zmf .* AREA_OCN;
    Amp_mean = Zmp .* AREA_OCN;
    Amd_mean = Zmd .* AREA_OCN;
    Alp_mean = Zlp .* AREA_OCN;
    Ald_mean = Zld .* AREA_OCN;
    Ab_mean  = Zb .* AREA_OCN;
    
    %% Calc LMEs
    for L=1:66
        lid = find(tlme==L);
        %total catch g
        lme_mcatch(L,t,1) = nansum(Amf_mcatch(lid));
        lme_mcatch(L,t,2) = nansum(Amp_mcatch(lid));
        lme_mcatch(L,t,3) = nansum(Amd_mcatch(lid));
        lme_mcatch(L,t,4) = nansum(Alp_mcatch(lid));
        lme_mcatch(L,t,5) = nansum(Ald_mcatch(lid));
        %mean biomass
        lme_mbio(L,t,1) = nanmean(Asf_mean(lid));
        lme_mbio(L,t,2) = nanmean(Asp_mean(lid));
        lme_mbio(L,t,3) = nanmean(Asd_mean(lid));
        lme_mbio(L,t,4) = nanmean(Amf_mean(lid));
        lme_mbio(L,t,5) = nanmean(Amp_mean(lid));
        lme_mbio(L,t,6) = nanmean(Amd_mean(lid));
        lme_mbio(L,t,7) = nanmean(Alp_mean(lid));
        lme_mbio(L,t,8) = nanmean(Ald_mean(lid));
        lme_mbio(L,t,9) = nanmean(Ab_mean(lid));
        %total biomass
        lme_sbio(L,t,1) = nansum(Asf_mean(lid));
        lme_sbio(L,t,2) = nansum(Asp_mean(lid));
        lme_sbio(L,t,3) = nansum(Asd_mean(lid));
        lme_sbio(L,t,4) = nansum(Amf_mean(lid));
        lme_sbio(L,t,5) = nansum(Amp_mean(lid));
        lme_sbio(L,t,6) = nansum(Amd_mean(lid));
        lme_sbio(L,t,7) = nansum(Alp_mean(lid));
        lme_sbio(L,t,8) = nansum(Ald_mean(lid));
        lme_sbio(L,t,9) = nansum(Ab_mean(lid));
        
    end
end

%%
save([fpath 'LME_core_',harv,'_' cfile '.mat'],...
    'lme_mcatch','lme_mbio','lme_sbio','lme_area');

%% Figures

%% Time series Plots 
% Global
figure(7)
subplot(2,2,1)
plot(years,sum(Flme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Flme_fcatch_all),'r','LineWidth',2)
ylabel('FEISTY catch (MT km^-^2)')
title('Forage')

subplot(2,2,2)
plot(years,sum(Plme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Plme_fcatch_all),'b','LineWidth',2)
title('Large pelagic')

subplot(2,2,3)
plot(years,sum(Dlme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Dlme_fcatch_all),'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
ylabel('FEISTY catch (MT km^-^2)')
title('Demersal')

subplot(2,2,4)
plot(years,sum(slme_scatch_all),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,sum(Alme_fcatch_all),'k','LineWidth',2)
xlabel('year')
title('All')
print('-dpng',[ppath 'CORE_',harv,'_sumlme_ts.png'])

%% log10 Global
figure(8)
subplot(2,2,1)
plot(years,log10(sum(Flme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Flme_fcatch_all)),'r','LineWidth',2)
ylabel('log_1_0 FEISTY catch (MT km^-^2)')
title('Forage')

subplot(2,2,2)
plot(years,log10(sum(Plme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Plme_fcatch_all)),'b','LineWidth',2)
title('Large pelagic')

subplot(2,2,3)
plot(years,log10(sum(Dlme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Dlme_fcatch_all)),'color',[0 0.6 0],'LineWidth',2)
xlabel('year')
ylabel('log_1_0 FEISTY catch (MT km^-^2)')
title('Demersal')

subplot(2,2,4)
plot(years,log10(sum(slme_scatch_all)),'color',[0.5 0.5 0.5],'LineWidth',2); hold on
plot(years,log10(sum(Alme_fcatch_all)),'k','LineWidth',2)
xlabel('year')
title('All')
print('-dpng',[ppath 'CORE_',harv,'_sumlme_ts_log10.png'])

%% TS by LME
% all 
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,slme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Alme_fcatch_all(i,:),'k','LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_',harv,'_ts_Allfish.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,Flme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Flme_fcatch_all(i,:),'r','LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_',harv,'_ts_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,Plme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Plme_fcatch_all(i,:),'b','LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_',harv,'_ts_P.png'])

%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(years,Dlme_scatch_all(i,:),'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        plot(years,Dlme_fcatch_all(i,:),'color',[0 0.6 0],'LineWidth',1)
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'CORE_byLME_',harv,'_ts_D.png'])


