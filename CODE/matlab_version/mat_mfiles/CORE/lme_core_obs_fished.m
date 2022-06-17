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

harv = 'fish_Fobs050_Pobs100_Dobs050'; %All_fish_obs
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

% clme_mf = NaN*ones(ni,nj);
% clme_mp = clme_mf;
% clme_md = clme_mf;
% clme_lp = clme_mf;
% clme_ld = clme_mf;
% 
% lme_sf = NaN*ones(ni,nj);
% lme_sp = lme_sf;
% lme_sd = lme_sf;
% lme_mf = lme_sf;
% lme_mp = lme_sf;
% lme_md = lme_sf;
% lme_lp = lme_sf;
% lme_ld = lme_sf;
% lme_b = lme_sf;
% 
% for L=1:66
%     lid = find(tlme==L);
%     
%     clme_mf(lid) = lme_mcatch(L,1);
%     clme_mp(lid) = lme_mcatch(L,2);
%     clme_md(lid) = lme_mcatch(L,3);
%     clme_lp(lid) = lme_mcatch(L,4);
%     clme_ld(lid) = lme_mcatch(L,5);
%     
%     lme_sf(lid) = lme_mbio(L,1);
%     lme_sp(lid) = lme_mbio(L,2);
%     lme_sd(lid) = lme_mbio(L,3);
%     lme_mf(lid) = lme_mbio(L,4);
%     lme_mp(lid) = lme_mbio(L,5);
%     lme_md(lid) = lme_mbio(L,6);
%     lme_lp(lid) = lme_mbio(L,7);
%     lme_ld(lid) = lme_mbio(L,8);
%     lme_b(lid) = lme_mbio(L,9);
% end
% 
% clme_All = clme_mf+clme_mp+clme_md+clme_lp+clme_ld;
% clme_AllF = clme_mf;
% clme_AllP = clme_mp+clme_lp;
% clme_AllD = clme_md+clme_ld;
% clme_AllM = clme_mf+clme_mp+clme_md;
% clme_AllL = clme_lp+clme_ld;
% 
% lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
% lme_AllF = lme_sf+lme_mf;
% lme_AllP = lme_sp+lme_mp+lme_lp;
% lme_AllD = lme_sd+lme_md+lme_ld;
% 
% %% Catch
% % all
% figure(1)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_All*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([5 7]);
% hcb = colorbar('h');
% set(gcf,'renderer','painters')
% title('Climatology LME mean log_1_0 total annual catch (MT) All Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_All.png'])
% 
% %% all F
% figure(2)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllF*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Forage Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllF.png'])
% 
% % all P
% figure(3)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllP*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Pelagic Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllP.png'])
% 
% % All D
% figure(4)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllD*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Demersal Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllD.png'])
% 
% % all M
% figure(5)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllM*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Medium Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllM.png'])
% 
% % all L
% figure(6)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% surfm(geolat_t,geolon_t,real(log10(clme_AllL*1e-6)))
% colormap('jet')
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% caxis([3.5 6.5]);
% hcb = colorbar('h');
% ylim(hcb,[3.5 6.5])                   %Set color axis if needed
% set(gcf,'renderer','painters')
% title('Historic 2051-2100 LME mean log10 total annual catch (MT) Large Fishes')
% stamp(cfile)
% print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch_AllL.png'])
% 
% 
% %% blank map
% figure(10)
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% print('-dpng',[ppath 'Map_blank.png'])

