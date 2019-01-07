% Estimate POEM harvest with different efforts scaled from land distance
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/CSV/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'esm26_1deg_shore_dist.mat']);

AREA_OCN = max(area,1);
tlme = lme_mask_onedeg;

cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_bio_prod_fish_Climatol_' harv '_' cfile '.mat']);

%% Plots in space
[ni,nj]=size(lon);

Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Cmf(ID)=mf_my;
Cmp(ID)=mp_my;
Cmd(ID)=md_my;
Clp(ID)=lp_my;
Cld(ID)=ld_my;

% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));



%% Effort scalings
mins = [0;0;0;0.5;0.5;0.5];
maxs = [1;1.5;2;1;1.5;2];

min_dist = -1*min_dist.^(1/3);
maxd = max(min_dist(:));
mind = min(min_dist(:));

%% apply effort scaling
for n=1:length(mins)
    
    eff = mins(n) + (maxs(n)-mins(n)).*(min_dist-mind) ./ (maxd-mind);
    % g/m2/d --> total g
    Amf_mcatch = Cmf .* AREA_OCN * 365 .* eff; %mean fish catch per yr
    Amp_mcatch = Cmp .* AREA_OCN * 365 .* eff;
    Amd_mcatch = Cmd .* AREA_OCN * 365 .* eff;
    Alp_mcatch = Clp .* AREA_OCN * 365 .* eff;
    Ald_mcatch = Cld .* AREA_OCN * 365 .* eff;
    
    %% Calc LMEs
    lme_mcatch = NaN*ones(66,5);
    
    for L=1:66
        lid = find(tlme==L);
        %total catch g
        lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
        lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
        lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
        lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
        lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
    end
    plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4));
    plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5));
    pFracPD = plme_Pmcatch ./ (plme_Pmcatch + plme_Dmcatch);
    
    %%
    cfile2 = ['_loop_effort_' num2str(mins(n)*10) num2str(maxs(n)*10)];
    save([fpath 'LME_clim_',harv,cfile2,'.mat'],'lme_mcatch');
    
    %% Figures
    
    clme_mf = NaN*ones(size(lon));
    clme_mp = clme_mf;
    clme_md = clme_mf;
    clme_lp = clme_mf;
    clme_ld = clme_mf;
    clme_PD = clme_mf;
    
    for L=1:66
        lid = find(tlme==L);
    
        clme_mf(lid) = lme_mcatch(L,1);
        clme_mp(lid) = lme_mcatch(L,2);
        clme_md(lid) = lme_mcatch(L,3);
        clme_lp(lid) = lme_mcatch(L,4);
        clme_ld(lid) = lme_mcatch(L,5);
        clme_PD(lid) = pFracPD(L);
    end
    
    clme_All = clme_mf+clme_mp+clme_md+clme_lp+clme_ld;
    clme_AllF = clme_mf;
    clme_AllP = clme_mp+clme_lp;
    clme_AllD = clme_md+clme_ld;
    clme_AllM = clme_mf+clme_mp+clme_md;
    clme_AllL = clme_lp+clme_ld;
    
    %% Catch
    
    % all
    figure(1)
    clf
    subplot('Position',[0.5 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,real(log10(clme_All*1e-6)))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([3 8]);
    set(gcf,'renderer','painters')
    title('All Fishes')
    
    % all F
    subplot('Position',[0 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,real(log10(clme_AllF*1e-6)))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([3 8]);
    %hcb = colorbar('h');
    colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
    set(gcf,'renderer','painters')
    title('mean log10 total annual catch (MT) Forage Fishes')
    
    % all P
    subplot('Position',[0.5 0.51 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,real(log10(clme_AllP*1e-6)))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([3 8]);
    %hcb = colorbar('h');
    set(gcf,'renderer','painters')
    title('Large Pelagic Fishes')
    
    % All D
    subplot('Position',[0 0 0.5 0.5])
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,real(log10(clme_AllD*1e-6)))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([3 8]);
    %hcb = colorbar('h');
    set(gcf,'renderer','painters')
    title('Demersal Fishes')
    stamp(cfile)
    print('-dpng',[ppath 'Clim_fished_',harv,'_LME_catch',cfile2,'.png'])
    
    
    figure(2)
    subplot(3,2,n)
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,clme_PD)
    cmocean('balance')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([0 1]);
    %hcb = colorbar('h');
    set(gcf,'renderer','painters')
    if (n==2)
        title('Fraction Large Pelagic vs. Demersal Catch')
    else
        title([num2str(mins(n)*10) '-' num2str(maxs(n)*10)])
    end
    
end
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_PDcatch_loop_effort.png'])
    

