% FEISTY F vs. RE pcolor of MMSSY
% Same F rate applied to all fish
% Climatology

clear all
close all

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE01000/'];

ppath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/',...
    'bio_rates/'];

load([dp 'Climatol_same_Frates_all_fish_RE_sens.mat']);

cmYOR=cbrewer('seq','YlOrRd',28);

%% total catches
lme_Amcatch = lme_Fmcatch + lme_Pmcatch + lme_Dmcatch;

Acatch = sum(MFc+MPc+MDc+LPc+LDc,3);
Fcatch = sum(MFc,3);
Pcatch = sum(MPc+LPc,3);
Dcatch = sum(MDc+LDc,3);

%% pcolor grids
jays = [Fish 7.4128];
kays = [reff 0.1668];
[kgrid,jgrid]=meshgrid(kays,jays);

nj = length(Fish);
nk = length(reff);

% Catch/Yield global & LME
cF2 = NaN*ones(nj+1,nk+1);
cP2 = NaN*ones(nj+1,nk+1);
cD2 = NaN*ones(nj+1,nk+1);
cA2 = NaN*ones(nj+1,nk+1);

lF2 = NaN*ones(nj+1,nk+1,66);
lP2 = NaN*ones(nj+1,nk+1,66);
lD2 = NaN*ones(nj+1,nk+1,66);
lA2 = NaN*ones(nj+1,nk+1,66);

cF2(1:nj,1:nk) = Fcatch;
cP2(1:nj,1:nk) = Pcatch;
cD2(1:nj,1:nk) = Dcatch;
cA2(1:nj,1:nk) = Acatch;

lF2(1:nj,1:nk,:) = lme_Fmcatch;
lP2(1:nj,1:nk,:) = lme_Pmcatch;
lD2(1:nj,1:nk,:) = lme_Dmcatch;
lA2(1:nj,1:nk,:) = lme_Amcatch;

% Comparison with SAUP
rF2 = NaN*ones(nj+1,nk+1);
rP2 = NaN*ones(nj+1,nk+1);
rD2 = NaN*ones(nj+1,nk+1);
rA2 = NaN*ones(nj+1,nk+1);

eF2 = NaN*ones(nj+1,nk+1);
eP2 = NaN*ones(nj+1,nk+1);
eD2 = NaN*ones(nj+1,nk+1);
eA2 = NaN*ones(nj+1,nk+1);

rF2(1:nj,1:nk) = r_all(:,:,2);
rP2(1:nj,1:nk) = r_all(:,:,3);
rD2(1:nj,1:nk) = r_all(:,:,4);
rA2(1:nj,1:nk) = r_all(:,:,1);

eF2(1:nj,1:nk) = rmse_all(:,:,2);
eP2(1:nj,1:nk) = rmse_all(:,:,3);
eD2(1:nj,1:nk) = rmse_all(:,:,4);
eA2(1:nj,1:nk) = rmse_all(:,:,1);


%% Plots of catch global
figure(1);
subplot(2,2,1)
pcolor(log(kgrid),log(jgrid),cF2)
colormap(cmYOR)
colorbar
%caxis([0 1])
ylabel('F (y^-^1)')
xlabel('RE')
title('F catch')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,2)
pcolor(log(kgrid),log(jgrid),cP2)
colormap(cmYOR)
colorbar
%caxis([0 1])
ylabel('F (y^-^1)')
xlabel('RE')
title('P catch')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,3)
pcolor(log(kgrid),log(jgrid),cD2)
colormap(cmYOR)
colorbar
%caxis([0 1])
ylabel('F (y^-^1)')
xlabel('RE')
title('D catch')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,4)
pcolor(log(kgrid),log(jgrid),cA2)
colormap(cmYOR)
colorbar
%caxis([0 1])
ylabel('F (y^-^1)')
xlabel('RE')
title('All catch')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))
print('-dpng',[ppath 'Climatol_global_yield_sameF_all_fish_RE_pcolor.png'])

%% Plots of catch by LME
nrows=11;
ncols=6;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);
ylab = 1:6:66;

% all
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        pcolor(log(kgrid),log(jgrid),lA2(:,:,i))
        colormap(cmYOR)
        %colorbar
        %caxis([0 1])
        if ~isempty(intersect(i,ylab))
            ylabel('F (y^-^1)')
            set(gca,'YTick',log(jays(2:3:end)),'YTickLabel',round(jays(1:3:end),2))
        end
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        else
            set(gca,'XTick',log(kays(2:4:end)),'XTickLabel',round(kays(1:4:end),3))
            xlabel('RE')
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_RE_pcolor_All.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        pcolor(log(kgrid),log(jgrid),lF2(:,:,i))
        colormap(cmYOR)
        %colorbar
        %caxis([0 1])
        if ~isempty(intersect(i,ylab))
            ylabel('F (y^-^1)')
            set(gca,'YTick',log(jays(2:3:end)),'YTickLabel',round(jays(1:3:end),2))
        end
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        else
            set(gca,'XTick',log(kays(2:4:end)),'XTickLabel',round(kays(1:4:end),3))
            xlabel('RE')
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_RE_pcolor_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        pcolor(log(kgrid),log(jgrid),lP2(:,:,i))
        colormap(cmYOR)
        %colorbar
        %caxis([0 1])
        if ~isempty(intersect(i,ylab))
            ylabel('F (y^-^1)')
            set(gca,'YTick',log(jays(2:3:end)),'YTickLabel',round(jays(1:3:end),2))
        end
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        else
            set(gca,'XTick',log(kays(2:4:end)),'XTickLabel',round(kays(1:4:end),3))
            xlabel('RE')
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_RE_pcolor_P.png'])


%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        pcolor(log(kgrid),log(jgrid),lD2(:,:,i))
        colormap(cmYOR)
        %colorbar
        %caxis([0 1])
        if ~isempty(intersect(i,ylab))
            ylabel('F (y^-^1)')
            set(gca,'YTick',log(jays(2:3:end)),'YTickLabel',round(jays(1:3:end),2))
        end
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        else
            set(gca,'XTick',log(kays(2:4:end)),'XTickLabel',round(kays(1:4:end),3))
            xlabel('RE')
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_RE_pcolor_D.png'])


%% Find MSY
fid = squeeze(max(max(lme_Fmcatch,[],1),[],2));
pid = squeeze(max(max(lme_Pmcatch,[],1),[],2));
did = squeeze(max(max(lme_Dmcatch,[],1),[],2));
aid = squeeze(max(max(lme_Amcatch,[],1),[],2));

[rgrid,fgrid]=meshgrid(reff,Fish);

FFmsy = NaN*ones(66,2);
PFmsy = NaN*ones(66,2);
DFmsy = NaN*ones(66,2);
AFmsy = NaN*ones(66,2);
for n=1:66
    fid2 = find(lme_Fmcatch(:,:,n)==fid(n));
    pid2 = find(lme_Pmcatch(:,:,n)==pid(n));
    did2 = find(lme_Dmcatch(:,:,n)==did(n));
    aid2 = find(lme_Amcatch(:,:,n)==aid(n));
    FFmsy(n,1) = rgrid(fid2);
    FFmsy(n,2) = fgrid(fid2);
    PFmsy(n,1) = rgrid(pid2);
    PFmsy(n,2) = fgrid(pid2);
    DFmsy(n,1) = rgrid(did2);
    DFmsy(n,2) = fgrid(did2);
    AFmsy(n,1) = rgrid(aid2);
    AFmsy(n,2) = fgrid(aid2);
end

%Fmsy5 = (Fmsy==5);
mean(FFmsy) %0.0596    3.4089
mean(PFmsy) %0.0739    4.0460
mean(DFmsy) %0.0055    4.8444
mean(AFmsy) %0.0663    4.8729

%% Plots of comp with SAUP
%Corr
figure(6);
subplot(2,2,1)
pcolor(log(kgrid),log(jgrid),rF2)
colormap(cmYOR)
colorbar
caxis([0.2 0.8])
ylabel('F (y^-^1)')
xlabel('RE')
title('F corr')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,2)
pcolor(log(kgrid),log(jgrid),rP2)
colormap(cmYOR)
colorbar
caxis([0.2 0.8])
ylabel('F (y^-^1)')
xlabel('RE')
title('P corr')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,3)
pcolor(log(kgrid),log(jgrid),rD2)
colormap(cmYOR)
colorbar
caxis([0.2 0.8])
ylabel('F (y^-^1)')
xlabel('RE')
title('D corr')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,4)
pcolor(log(kgrid),log(jgrid),rA2)
colormap(cmYOR)
colorbar
caxis([0.2 0.8])
ylabel('F (y^-^1)')
xlabel('RE')
title('All corr')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))
print('-dpng',[ppath 'Climatol_global_SAUPcorr_sameF_all_fish_RE_pcolor.png'])

%% RMSE
cmROY = flipud(cmYOR);

figure(7);
subplot(2,2,1)
pcolor(log(kgrid),log(jgrid),eF2)
colormap(cmROY)
colorbar
%caxis([0 2])
ylabel('F (y^-^1)')
xlabel('RE')
title('F RMSE')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,2)
pcolor(log(kgrid),log(jgrid),eP2)
colormap(cmROY)
colorbar
%caxis([0 2])
ylabel('F (y^-^1)')
xlabel('RE')
title('P RMSE')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,3)
pcolor(log(kgrid),log(jgrid),eD2)
colormap(cmROY)
colorbar
%caxis([0 2])
ylabel('F (y^-^1)')
xlabel('RE')
title('D RMSE')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))

subplot(2,2,4)
pcolor(log(kgrid),log(jgrid),eA2)
colormap(cmROY)
colorbar
%caxis([0 2])
ylabel('F (y^-^1)')
xlabel('RE')
title('All RMSE')
set(gca,'XTick',log(kays(2:2:end)),'XTickLabel',round(kays(1:2:end),3),...
    'YTick',log(Fish(2:2:end)),'YTickLabel',round(jays(1:2:end),2))
print('-dpng',[ppath 'Climatol_global_SAUPrmse_sameF_all_fish_RE_pcolor.png'])

%% Find max corr & min RMSE
fid = squeeze(max(max(r_all(:,:,2),[],1),[],2));
pid = squeeze(max(max(r_all(:,:,3),[],1),[],2));
did = squeeze(max(max(r_all(:,:,4),[],1),[],2));
aid = squeeze(max(max(r_all(:,:,1),[],1),[],2));
fid2 = find(r_all(:,:,2)==fid);
pid2 = find(r_all(:,:,3)==pid);
did2 = find(r_all(:,:,4)==did);
aid2 = find(r_all(:,:,1)==aid);
maxR(1,1) = rgrid(fid2);
maxR(1,2) = fgrid(fid2);
maxR(2,1) = rgrid(pid2);
maxR(2,2) = fgrid(pid2);
maxR(3,1) = rgrid(did2);
maxR(3,2) = fgrid(did2);
maxR(4,1) = rgrid(aid2);
maxR(4,2) = fgrid(aid2);

fie = squeeze(min(min(rmse_all(:,:,2),[],1),[],2));
pie = squeeze(min(min(rmse_all(:,:,3),[],1),[],2));
die = squeeze(min(min(rmse_all(:,:,4),[],1),[],2));
aie = squeeze(min(min(rmse_all(:,:,1),[],1),[],2));
fie2 = find(rmse_all(:,:,2)==fie);
pie2 = find(rmse_all(:,:,3)==pie);
die2 = find(rmse_all(:,:,4)==die);
aie2 = find(rmse_all(:,:,1)==aie);
minE(1,1) = rgrid(fie2);
minE(1,2) = fgrid(fie2);
minE(2,1) = rgrid(pie2);
minE(2,2) = fgrid(pie2);
minE(3,1) = rgrid(die2);
minE(3,2) = fgrid(die2);
minE(4,1) = rgrid(aie2);
minE(4,2) = fgrid(aie2);

