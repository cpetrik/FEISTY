% Visualize difference between
% CESM-WACCM Historic 1949 Spinup w/fishing and
% CESM FOSI 1948 Spinup w/fishing

% MAYBE NEED TO AGGREGATE BY REGIONS?

clear
close all

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
HS = Hsp+Hsf+Hsd;
HM = Hmp+Hmf+Hmd;
HL = Hlp+Hld;

%% Interpolate to same grid

hF = griddata(tlat,tlon,HF,glat,glon);
hP = griddata(tlat,tlon,HP,glat,glon);
hD = griddata(tlat,tlon,HD,glat,glon);
hB = griddata(tlat,tlon,Hb,glat,glon);
hS = griddata(tlat,tlon,HS,glat,glon);
hM = griddata(tlat,tlon,HM,glat,glon);
hL = griddata(tlat,tlon,HL,glat,glon);

hAll = hF+hP+hD;
hFrahPD = hP ./ (hP+hD);
hFrahPF = hP ./ (hP+hF);
hFrahLM = hL ./ (hL+hM);

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
    CS = Csp+Csf+Csd;
    CM = Cmp+Cmf+Cmd;
    CL = Clp+Cld;

    %% Interpolate to same grid

    cF = griddata(CLAT,CLON,CF,glat,glon);
    cP = griddata(CLAT,CLON,CP,glat,glon);
    cD = griddata(CLAT,CLON,CD,glat,glon);
    cB = griddata(CLAT,CLON,Cb,glat,glon);
    cS = griddata(CLAT,CLON,CS,glat,glon);
    cM = griddata(CLAT,CLON,CM,glat,glon);
    cL = griddata(CLAT,CLON,CL,glat,glon);

    %
    cAll = cF+cP+cD;
    cFracPD = cP ./ (cP+cD);
    cFracPF = cP ./ (cP+cF);
    cFracLM = cL ./ (cL+cM);

    %%
    diffF = (cF - hF);
    diffP = (cP - hP);
    diffD = (cD - hD);
    diffB = (cB - hB);
    diffAll = (cAll - hAll);

    pdiffF = (cF-hF) ./ hF;
    pdiffP = (cP-hP) ./ hP;
    pdiffD = (cD-hD) ./ hD;
    pdiffB = (cB-hB) ./ hB;
    pdiffAll = (cAll-hAll) ./ hAll;

    %% find non-nan grid cells (ocean cells)
    cid = find(~isnan(cAll));
    hid = find(~isnan(hAll));
    keep = intersect(cid,hid);

    %% Stats
    %Pearson correlation coeff (r)
    [rall,pall] =corr(cAll(keep),hAll(keep));
    [rF,pF] =corr(cF(keep),hF(keep));
    [rP,pP] =corr(cP(keep),hP(keep));
    [rD,pD] =corr(cD(keep),hD(keep));
    [rPD,pPD] =corr(cFracPD(keep),hFrahPD(keep));


    %Concordance correlation coeff (CCC)
    aCCC = f_CCC([cAll(keep),hAll(keep)],0.05);
    fCCC = f_CCC([cF(keep),hF(keep)],0.05);
    pCCC = f_CCC([cP(keep),hP(keep)],0.05);
    dCCC = f_CCC([cD(keep),hD(keep)],0.05);
    pdCCC = f_CCC([cFracPD(keep),hFrahPD(keep)],0.05);


    %root mean square error
    o=hAll(keep);
    p=cAll(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmse = sqrt(num/n);

    o=hF(keep);
    p=cF(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmseF = sqrt(num/n);

    o=hP(keep);
    p=cP(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmseP = sqrt(num/n);

    o=hD(keep);
    p=cD(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmseD = sqrt(num/n);

    o=hFrahPD(keep);
    p=cFracPD(keep);
    n = length(o);
    num=nansum((p-o).^2);
    rmsePD = sqrt(num/n);

    %Fmed
    Fall=10^(median(hAll(keep)-cAll(keep)));
    FF=10^(median(hF(keep)-cF(keep)));
    FP=10^(median(hP(keep)-cP(keep)));
    FD=10^(median(hD(keep)-cD(keep)));
    FPD=10^(median(hFrahPD(keep)-cFracPD(keep)));

    % Bias (FOSI minus SAUP)
    %average error = bias
    p=hAll(keep);
    o=cAll(keep);
    n = length(o);
    bias = nansum(o-p) / n;

    p=hF(keep);
    o=cF(keep);
    n = length(o);
    biasF = nansum(o-p) / n;

    p=hP(keep);
    o=cP(keep);
    n = length(o);
    biasP = nansum(o-p) / n;

    p=hD(keep);
    o=cD(keep);
    n = length(o);
    biasD = nansum(o-p) / n;

    p=hFrahPD(keep);
    o=cFracPD(keep);
    n = length(o);
    biasPD = nansum(o-p) / n;

    
    %MEF
    %(sum((r - rbar).^2) - sum((bsxfun(@minus, [r f], r)).^2))./(sum((r - rbar).^2));
    o=hAll(keep);
    p=cAll(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mef = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hF(keep);
    p=cF(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefF = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hP(keep);
    p=cP(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefP = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hD(keep);
    p=cD(keep);
    obar = nanmean(o);
    pbar = nanmean(p);
    mefD = (sum((o - obar).^2) - sum((p - pbar).^2)) ./ (sum((o - obar).^2));

    o=hFrahPD(keep);
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
    %scatter(StockPNAS(:,7),glme_mcatch(keep),20,lme_tp_fosi(keep,1),'filled'); hold on;
    %cmocean('thermal');
    scatter(log10(hF(keep)),log10(cF(keep)),20,'filled'); hold on;
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
    scatter(log10(hP(keep)),log10(cP(keep)),20,'filled'); hold on;
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
    plot(x,x5h,':r'); hold on;
    plot(x,x5l,':r'); hold on;
    scatter((hFrahPD(keep)),(cFracPD(keep)),20,'filled'); hold on;
    text(-0.05,1.0,['r = ' sprintf('%2.2f',rPD) ' '])
    text(-0.05,0.75,['RMSE = ' sprintf('%2.2f',rmsePD)])
    text(-0.05,0.5,['bias = ' sprintf('%2.2f',biasPD)])
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
    scatter(log10(hAll(keep)),log10(cAll(keep)),20,'filled'); hold on;
    text(-0.75,1.75,['r = ' sprintf('%2.2f',rall) ' '])
    text(-0.75,1.5,['RMSE = ' sprintf('%2.2f',rmse)])
    text(-0.75,1.25,['bias = ' sprintf('%2.2f',bias)])
    axis([-1 2.25 -1 2.25])
    xlabel('CESM FOSI')
    ylabel('CESM-WACCM')
    title('All fishes')
    stamp('')
    print('-dpng',[ppath 'WACCM_FOSI_' mod 'scatter_types.png'])

    %% Create a table for stats

end

