% Calc Fish-MIP outputs and save as NetCDF
% Preindustrial simulation

clear all
close all

%% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/FishMIP_CESM/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'pristine';
tharv = 'F=0';
fpath=['/Volumes/GFDL/NC/FishMIP/CESM1-BEC/' cfile '/'];

load([fpath 'Means_Preindust_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean',...
    'mf_mean','mp_mean','md_mean',...
    'lp_mean','ld_mean','b_mean');

ppath=['/Volumes/GFDL/Fish-MIP/CESM/PreIndust/'];
load([ppath 'cesm_pi_phy_zint_monthly_185001-210012.mat'])
load([ppath 'cesm_pi_zoo_zint_monthly_185001-210012.mat'])

Years = 1850:2100;

AllF = sf_mean + mf_mean;
AllP = sp_mean + mp_mean + lp_mean;
AllD = sd_mean + md_mean + ld_mean;
AllM = mp_mean + mf_mean + md_mean;
AllL = lp_mean + ld_mean;
All  = AllF + AllP + AllD;

%% Reshape to ocean cells only
cpath = '/Volumes/GFDL/Fish-MIP/CESM/';
load([cpath 'gridspec_cesm.mat'],'LAT','LON');
load([cpath 'Data_grid_cesm.mat']);

[ni,nj,NT] = size(nz_int);

phy = reshape(np_int,ni*nj,NT);
zoo = reshape(nz_int,ni*nj,NT);

vphy = phy(GRD.ID,:);
vzoo = zoo(GRD.ID,:);

clear np_int nz_int phy zoo

%% Units
%fish: gWW/m^2
%zoo_int: mmol C/m^2
%phy_int: mmol C/m^2
%tp: degC
%tb: degC

% Convert to gC m-2 from gWW m-2
AllM = (1/9) * AllM;
AllL = (1/9) * AllL;
All  = (1/9) * All;

% Convert to gC m-2 from mmol C m-2
% 1e-3 mol in 1 mmol
% 12.01 g C in 1 mol C
zint = vphy * 1e-3 * 12.01;
pint = vzoo * 1e-3 * 12.01;

%% Annual means of phy and zoo
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
t=1:NT;
mo=t/12;
mo=mo+1850;
st=1:12:NT;
en=12:12:NT;

phy_mean = nan(size(sd_mean));
zoo_mean = nan(size(sd_mean));
for m=1:length(en)
    yr1 = st(m):en(m);
    phy_mean(:,m)=mean(vphy(:,yr1),2);
    zoo_mean(:,m)=mean(vzoo(:,yr1),2);
end

%% Outputs
%Total system carbon biomass, tsb, gCm-2, All primary producers and consumers
%Total consumer carbon biomass density, tcb, gCm-2 All consumers (trophic level >1, vertebrates and invertebrates)
%Carbon biomass density of consumers > 10cm, b10cm, gCm-2 If asymptotic length (Linf) is > 10cm, include in > 10cm class
%Carbon biomass density of consumers > 30cm, b30cm, gCm-2 If asymptotic length (Linf) is > 30cm, include in > 30cm class
%Monthly or annual

% Ryan H said to change to only outputs from my model, so no phyto or zoop
%tsb = phy + zoo + F + P + D + B
%tcb = zoo + F + P + D + B
%b10cm = M + L + (0.1-0.25)*B
%b30cm = L

vtsb = All + b_mean;
vtcb = All + b_mean;
vb10cm = AllM + AllL + (0.1)*b_mean;
vb30cm = AllL;

%% Reshape to lat,lon,yr
[nid,nt] = size(sd_mean);

tsb = 1.0e36*ones(ni,nj,nt);
tcb = 1.0e36*ones(ni,nj,nt);
b10cm = 1.0e36*ones(ni,nj,nt);
b30cm = 1.0e36*ones(ni,nj,nt);

for y=1:length(en)
    gtsb = nan(ni,nj);
    ttsb = vtsb(:,y);
    gtsb(GRD.ID) = ttsb;
    tsb(:,:,y) = gtsb;
    
    gtcb = nan(ni,nj);
    ttcb = vtcb(:,y);
    gtcb(GRD.ID) = ttcb;
    tcb(:,:,y) = gtcb;
    
    gb10cm = nan(ni,nj);
    tb10cm = vb10cm(:,y);
    gb10cm(GRD.ID) = tb10cm;
    b10cm(:,:,y) = gb10cm;
    
    gb30cm = nan(ni,nj);
    tb30cm = vb30cm(:,y);
    gb30cm(GRD.ID) = tb30cm;
    b30cm(:,:,y) = gb30cm;
end

save([fpath 'FishMIP_output_Preindust_' cfile '.mat'],'Years',...
    'tsb','tcb','b10cm','b30cm','LAT','LON');

%% Quick look
sb = tsb(:,:,150);
cb = tcb(:,:,150);
b10 = b10cm(:,:,150);
b30 = b30cm(:,:,150);

figure(1)
pcolor(log10(sb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])

figure(2)
pcolor(log10(cb'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])

figure(3)
pcolor(log10(b10'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])

figure(4)
pcolor(log10(b30'))
shading flat
colormap('jet')
colorbar
caxis([-2 2])