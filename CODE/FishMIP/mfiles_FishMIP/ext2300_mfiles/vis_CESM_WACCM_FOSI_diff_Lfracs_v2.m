% Visualize difference between
% CESM-WACCM Historic 1949 Spinup w/fishing and
% CESM FOSI 1948 Spinup w/fishing

% MAYBE NEED TO AGGREGATE BY REGIONS?
%try using GFDL 1 degree LME map from ISIMIP
%with obsGLMM biomes for high seas

clear
close all

%% LME mask
load(['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/',...
    'lme_gfdl-mom6-cobalt2_onedeg.mat'], 'tlme')

figure
pcolor(tlme); shading flat

lme_mask = fliplr(tlme);

%% Biome mask
load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/biome_masks/ESM_Biome_Masks/data_biomes_MODISAqua_x1.mat',...
    'biomes')

figure
pcolor(biomes); shading flat

%% CESM-WACCM
wpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/CESM2-WACCM/';
load([wpath 'gridspec_cesm2_cmip6_2300.mat']);
load([wpath 'Data_grid_cesm2_cmip6_2300.mat']);

CID = GRD.ID;
CLAT = LAT;
CLON = LON;

clear GRD LAT LON

[hi,hj]=size(CLON);

%% CESM FOSI grid
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat'],'TLAT','TLONG');
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);
FID = GRD.ID;

%% Need to fix FOSI longitude
%tlat   [-79.2205 89.7064]
%clat   [-89.5 89.5]
%tlon   [0.0147 359.996]
%clon   [-179.5 179.5]

test = TLONG;
id=find(test>180);
test(id)=test(id)-360;
tlon = test;
tlat = TLAT;

lats = -89.5:89.5;
lons = -179.5:179.5;
[glon,glat] = meshgrid(lons,lats);

glon=glon';
glat=glat';

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

load coastlines

%% FEISTY Output
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
cfile2 = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
harv = 'All_fish03';

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%% Put FOSI results outside loop
% CESM-FOSI ======================================================
% /Volumes/petrik-lab/Feisty/NC/CESM_MAPP/cfile/Spinup/
%Means_Spinup_v15_ cfile
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile2 '/Spinup/'];
load([fpath 'Means_Spinup_v15_' cfile2 '.mat']);

%% Put biomass on grid
Hsf=NaN*ones(ni,nj);
Hsp=NaN*ones(ni,nj);
Hsd=NaN*ones(ni,nj);
Hmf=NaN*ones(ni,nj);
Hmp=NaN*ones(ni,nj);
Hmd=NaN*ones(ni,nj);
Hlp=NaN*ones(ni,nj);
Hld=NaN*ones(ni,nj);
Hb =NaN*ones(ni,nj);

Hsf(FID)=sf_mean;
Hsp(FID)=sp_mean;
Hsd(FID)=sd_mean;
Hmf(FID)=mf_mean;
Hmp(FID)=mp_mean;
Hmd(FID)=md_mean;
Hlp(FID)=lp_mean;
Hld(FID)=ld_mean;
Hb(FID) =b_mean;

Hsf_ts=sf_tmean;
Hsp_ts=sp_tmean;
Hsd_ts=sd_tmean;
Hmf_ts=mf_tmean;
Hmp_ts=mp_tmean;
Hmd_ts=md_tmean;
Hlp_ts=lp_tmean;
Hld_ts=ld_tmean;
Hb_ts =b_tmean;

%%
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean lp_tmean ld_tmean b_tmean
clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean b_mean

HF = Hsf+Hmf;
HP = Hsp+Hmp+Hlp;
HD = Hsd+Hmd+Hld;

%% Interpolate to same grid

hF = griddata(tlat,tlon,HF,glat,glon);
hP = griddata(tlat,tlon,HP,glat,glon);
hD = griddata(tlat,tlon,HD,glat,glon);

hAll = hF+hP+hD;

%% CESM-WACCM ======================================================
gpath=['/Volumes/petrik-lab/Feisty/NC/WG2300/',cfile,'/CESM2-WACCM/'];
fracs = {'50','55','60','65','70','75','80','85','90','95','zooc'};

%%
for m = 9:length(fracs)
    mod = [fracs{m} '_'];
    if m==length(fracs)
        exper = ['CESM2-WACCM_spinup_',mod,'All_fish03'];
    else
        exper = ['CESM2-WACCM_spinup_Lfrac',mod,'All_fish03'];
    end

    load([gpath 'Means_' exper '_' cfile '.mat']);

    close all

    %%
    Csf=NaN*ones(hi,hj);
    Csp=NaN*ones(hi,hj);
    Csd=NaN*ones(hi,hj);
    Cmf=NaN*ones(hi,hj);
    Cmp=NaN*ones(hi,hj);
    Cmd=NaN*ones(hi,hj);
    Clp=NaN*ones(hi,hj);
    Cld=NaN*ones(hi,hj);
    Cb =NaN*ones(hi,hj);

    Csf(CID)=sf_mean;
    Csp(CID)=sp_mean;
    Csd(CID)=sd_mean;
    Cmf(CID)=mf_mean;
    Cmp(CID)=mp_mean;
    Cmd(CID)=md_mean;
    Clp(CID)=lp_mean;
    Cld(CID)=ld_mean;
    Cb(CID) =b_mean;

    Csf_ts=sf_tmean;
    Csp_ts=sp_tmean;
    Csd_ts=sd_tmean;
    Cmf_ts=mf_tmean;
    Cmp_ts=mp_tmean;
    Cmd_ts=md_tmean;
    Clp_ts=lp_tmean;
    Cld_ts=ld_tmean;
    Cb_ts =b_tmean;

    %%
    clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean lp_tmean ld_tmean b_tmean
    clear sf_mean sp_mean sd_mean mf_mean mp_mean md_mean lp_mean ld_mean b_mean

    %%
    CF = Csf+Cmf;
    CP = Csp+Cmp+Clp;
    CD = Csd+Cmd+Cld;
    
    %% Interpolate to same grid

    cF = griddata(CLAT,CLON,CF,glat,glon);
    cP = griddata(CLAT,CLON,CP,glat,glon);
    cD = griddata(CLAT,CLON,CD,glat,glon);
    
    %
    cAll = cF+cP+cD;
    
    %% LME means
    % should be area-weighted means but ignore for now

    clmeAll = NaN*ones(69,1);
    clmeF = NaN*ones(69,1);
    clmeP = NaN*ones(69,1);
    clmeD = NaN*ones(69,1);

    hlmeAll = NaN*ones(69,1);
    hlmeF = NaN*ones(69,1);
    hlmeP = NaN*ones(69,1);
    hlmeD = NaN*ones(69,1);

    for L=1:66
        lid = find(lme_mask==L);
        % waccm
        clmeAll(L,1) = mean(cAll(lid),'omitnan');
        clmeF(L,1) = mean(cF(lid),'omitnan');
        clmeP(L,1) = mean(cP(lid),'omitnan');
        clmeD(L,1) = mean(cD(lid),'omitnan');

        % fosi
        hlmeAll(L,1) = mean(hAll(lid),'omitnan');
        hlmeF(L,1) = mean(hF(lid),'omitnan');
        hlmeP(L,1) = mean(hP(lid),'omitnan');
        hlmeD(L,1) = mean(hD(lid),'omitnan');
    end

    %find biomes in high seas
    nid = find(isnan(lme_mask));
    for b=1:3
        hid = find(biomes==b);
        bid = intersect(hid,nid);

        % waccm
        clmeAll(66+b,1) = mean(cAll(bid),'omitnan');
        clmeF(66+b,1) = mean(cF(bid),'omitnan');
        clmeP(66+b,1) = mean(cP(bid),'omitnan');
        clmeD(66+b,1) = mean(cD(bid),'omitnan');

        % fosi
        hlmeAll(66+b,1) = mean(hAll(bid),'omitnan');
        hlmeF(66+b,1) = mean(hF(bid),'omitnan');
        hlmeP(66+b,1) = mean(hP(bid),'omitnan');
        hlmeD(66+b,1) = mean(hD(bid),'omitnan');
    end

    %% PD frac
    hFracPD = hlmeP ./ (hlmeP+hlmeD);
    cFracPD = clmeP ./ (clmeP+clmeD);
    
    %remove interior seas 
    keep = [1:22,24:32,34:61,63:69];

    %% Stats
    %Pearson correlation coeff (r)
    [rall,pall] =corr(clmeAll(keep),hlmeAll(keep));
    [rF,pF] =corr(clmeF(keep),hlmeF(keep));
    [rP,pP] =corr(clmeP(keep),hlmeP(keep));
    [rD,pD] =corr(clmeD(keep),hlmeD(keep));
    [rPD,pPD] =corr(cFracPD(keep),hFracPD(keep));


    %Concordance correlation coeff (CCC)
    aCCC = f_CCC([clmeAll(keep),hlmeAll(keep)],0.05);
    fCCC = f_CCC([clmeF(keep),hlmeF(keep)],0.05);
    pCCC = f_CCC([clmeP(keep),hlmeP(keep)],0.05);
    dCCC = f_CCC([clmeD(keep),hlmeD(keep)],0.05);
    pdCCC = f_CCC([cFracPD(keep),hFracPD(keep)],0.05);


    %root mean square error
    o=hlmeAll(keep);
    p=clmeAll(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse = sqrt(num/n);

    o=hlmeF(keep);
    p=clmeF(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmseF = sqrt(num/n);

    o=hlmeP(keep);
    p=clmeP(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmseP = sqrt(num/n);

    o=hlmeD(keep);
    p=clmeD(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmseD = sqrt(num/n);

    o=hFracPD(keep);
    p=cFracPD(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmsePD = sqrt(num/n);

    %Fmed
    Fall=10^(median(hlmeAll(keep)-clmeAll(keep)));
    FF=10^(median(hlmeF(keep)-clmeF(keep)));
    FP=10^(median(hlmeP(keep)-clmeP(keep)));
    FD=10^(median(hlmeD(keep)-clmeD(keep)));
    FPD=10^(median(hFracPD(keep)-cFracPD(keep)));

    % Bias (FOSI minus SAUP)
    %average error = bias
    p=hlmeAll(keep);
    o=clmeAll(keep);
    n = length(o);
    bias = nansum(o-p) / n;

    p=hlmeF(keep);
    o=clmeF(keep);
    n = length(o);
    biasF = nansum(o-p) / n;

    p=hlmeP(keep);
    o=clmeP(keep);
    n = length(o);
    biasP = nansum(o-p) / n;

    p=hlmeD(keep);
    o=clmeD(keep);
    n = length(o);
    biasD = nansum(o-p) / n;

    p=hFracPD(keep);
    o=cFracPD(keep);
    n = length(o);
    biasPD = nansum(o-p) / n;


    %MEF
    %(sum((r - rbar).^2) - sum((bsxfun(@minus, [r f], r)).^2))./(sum((r - rbar).^2));
    o=hlmeAll(keep);
    p=clmeAll(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mef = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hlmeF(keep);
    p=clmeF(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefF = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hlmeP(keep);
    p=clmeP(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefP = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hlmeD(keep);
    p=clmeD(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefD = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hFracPD(keep);
    p=cFracPD(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefPD = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    %% Scatter plots

    x=-5:0.5:5;
    x2h = x+log10(2);
    x2l = x-log10(2);
    x5h = x+log10(5);
    x5l = x-log10(5);

    figure(1)
    subplot(2,2,1)
    plot(x,x,'--k'); hold on;
    plot(x,x2h,':b'); hold on;
    plot(x,x2l,':b'); hold on;
    plot(x,x5h,':r'); hold on;
    plot(x,x5l,':r'); hold on;
    scatter(log10(hlmeF(keep)),log10(clmeF(keep)),20,'filled'); hold on;
    text(-3.25,1.25,['r = ' sprintf('%2.2f',rF) ' '])
    text(-3.25,0.75,['RMSE = ' sprintf('%2.2f',rmseF)])
    text(-3.25,0.25,['bias = ' sprintf('%2.2f',biasF)])
    axis([-4 1.5 -4 1.5])
    xlabel('CESM FOSI')
    ylabel('CESM-WACCM')
    title('Forage fishes')

    subplot(2,2,2)
    plot(x,x,'--k'); hold on;
    plot(x,x2h,':b'); hold on;
    plot(x,x2l,':b'); hold on;
    plot(x,x5h,':r'); hold on;
    plot(x,x5l,':r'); hold on;
    scatter(log10(hlmeP(keep)),log10(clmeP(keep)),20,'filled'); hold on;
    text(-3.25,1.25,['r = ' sprintf('%2.2f',rP) ' '])
    text(-3.25,0.75,['RMSE = ' sprintf('%2.2f',rmseP)])
    text(-3.25,0.25,['bias = ' sprintf('%2.2f',biasP)])
    axis([-4 1.5 -4 1.5])
    xlabel('CESM FOSI')
    ylabel('CESM-WACCM')
    title('LgPel fishes')

    subplot(2,2,3)
    plot(x,x,'--k'); hold on;
    plot(x,x2h,':b'); hold on;
    plot(x,x2l,':b'); hold on;
    scatter((hFracPD(keep)),(cFracPD(keep)),20,'filled'); hold on;
    text(-0.05,1.0,['r = ' sprintf('%2.2f',rPD) ' '])
    text(-0.05,0.9,['RMSE = ' sprintf('%2.2f',rmsePD)])
    text(-0.05,0.8,['bias = ' sprintf('%2.2f',biasPD)])
    axis([-0.1 1.1 -0.1 1.1])
    xlabel('CESM FOSI')
    ylabel('CESM-WACCM')
    title('P / (P+D)')

    subplot(2,2,4)
    plot(x,x,'--k'); hold on;
    plot(x,x2h,':b'); hold on;
    plot(x,x2l,':b'); hold on;
    plot(x,x5h,':r'); hold on;
    plot(x,x5l,':r'); hold on;
    scatter(log10(hlmeAll(keep)),log10(clmeAll(keep)),20,'filled'); hold on;
    text(-0.75,1.75,['r = ' sprintf('%2.2f',rall) ' '])
    text(-0.75,1.5,['RMSE = ' sprintf('%2.2f',rmse)])
    text(-0.75,1.25,['bias = ' sprintf('%2.2f',bias)])
    axis([-1 2.25 -1 2.25])
    xlabel('CESM FOSI')
    ylabel('CESM-WACCM')
    title('All fishes')
    stamp(fracs{m})
    print('-dpng',[ppath 'WACCM_FOSI_' mod 'scatter_types.png'])

    
end

