% FEISTY assim & aMet search 
% pcolor of corr & RMSE w/SAUP

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/bio_rates/';
ppath = [pp 'Lambda_amet_search_1meso_Dc_enc70-b200_m-b175-k086_c20-b250_D080_A050_nmort1_BE10_noCC_RE00100/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

dp = '/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/';
load([dp 'Lambda_amet_search_1meso_Dc_enc70-b200_m-b175-k086_c20-b250_D080_A050_nmort1_BE10_noCC_RE00100.mat'])

lam = 0.55:0.025:0.675;
aM = 4.25:0.25:6;

cmYOR=cbrewer('seq','YlOrRd',28,'PCHIP');

%% pcolor grids
lams = [lam 0.7];
mets = [aM 6.25];
[lgrid,mgrid]=meshgrid(lams,mets);

nj = length(lam);
nk = length(aM);

% Reshape
r_all = reshape(r_all,nk,nj,6);
rmse_all = reshape(rmse_all,nk,nj,6);

% Comparison with SAUP
rF2 = NaN*ones(nk+1,nj+1);
rP2 = NaN*ones(nk+1,nj+1);
rD2 = NaN*ones(nk+1,nj+1);
rA2 = NaN*ones(nk+1,nj+1);
rRS2 = NaN*ones(nk+1,nj+1);
rRD2 = NaN*ones(nk+1,nj+1);

eF2 = NaN*ones(nk+1,nj+1);
eP2 = NaN*ones(nk+1,nj+1);
eD2 = NaN*ones(nk+1,nj+1);
eA2 = NaN*ones(nk+1,nj+1);
eRS2 = NaN*ones(nk+1,nj+1);
eRD2 = NaN*ones(nk+1,nj+1);

rF2(1:nk,1:nj) = r_all(:,:,2);
rP2(1:nk,1:nj) = r_all(:,:,3);
rD2(1:nk,1:nj) = r_all(:,:,4);
rA2(1:nk,1:nj) = r_all(:,:,1);
rRS2(1:nk,1:nj) = r_all(:,:,5);
rRD2(1:nk,1:nj) = r_all(:,:,6);

eF2(1:nk,1:nj) = rmse_all(:,:,2);
eP2(1:nk,1:nj) = rmse_all(:,:,3);
eD2(1:nk,1:nj) = rmse_all(:,:,4);
eA2(1:nk,1:nj) = rmse_all(:,:,1);
eRS2(1:nk,1:nj) = rmse_all(:,:,5);
eRD2(1:nk,1:nj) = rmse_all(:,:,6);


%% Plots of comp with SAUP
%Corr
figure(1);
subplot(2,3,1)
pcolor(mgrid,lgrid,rF2)
colormap(cmYOR)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('F corr')

%
subplot(2,3,2)
pcolor(mgrid,lgrid,rP2)
colormap(cmYOR)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('P corr')


subplot(2,3,3)
pcolor(mgrid,lgrid,rD2)
colormap(cmYOR)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('D corr')


subplot(2,3,4)
pcolor(mgrid,lgrid,rD2)
colormap(cmYOR)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('All corr')

subplot(2,3,5)
pcolor(mgrid,lgrid,rRS2)
colormap(cmYOR)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('P/(P+D) SAU corr')

subplot(2,3,6)
pcolor(mgrid,lgrid,rRD2)
colormap(cmYOR)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('P/(P+D) DvD corr')

print('-dpng',[ppath 'Hist_global_SAUPcorr_assim_amet_pcolor.png'])

%% RMSE
cmROY = flipud(cmYOR);

figure(2);
subplot(2,3,1)
pcolor(mgrid,lgrid,eF2)
colormap(cmROY)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('F RMSE')

subplot(2,3,2)
pcolor(mgrid,lgrid,eP2)
colormap(cmROY)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('P RMSE')

subplot(2,3,3)
pcolor(mgrid,lgrid,eD2)
colormap(cmROY)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('D RMSE')

subplot(2,3,4)
pcolor(mgrid,lgrid,eA2)
colormap(cmROY)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('All RMSE')

subplot(2,3,5)
pcolor(mgrid,lgrid,eRS2)
colormap(cmROY)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('P/(P+D) SAUP corr')

subplot(2,3,6)
pcolor(mgrid,lgrid,eRD2)
colormap(cmROY)
colorbar
axis square
ylabel('assim eff')
xlabel('a met')
title('P/(P+D) DvD corr')

print('-dpng',[ppath 'Hist_global_SAUPrmse_assim_amet_pcolor.png'])

%% Find max corr & min RMSE
fid = squeeze(max(max(r_all(:,:,2),[],1),[],2));
pid = squeeze(max(max(r_all(:,:,3),[],1),[],2));
did = squeeze(max(max(r_all(:,:,4),[],1),[],2));
aid = squeeze(max(max(r_all(:,:,1),[],1),[],2));
sid = squeeze(max(max(r_all(:,:,5),[],1),[],2));
rid = squeeze(max(max(r_all(:,:,6),[],1),[],2));
fid2 = find(r_all(:,:,2)==fid);
pid2 = find(r_all(:,:,3)==pid);
did2 = find(r_all(:,:,4)==did);
aid2 = find(r_all(:,:,1)==aid);
sid2 = find(r_all(:,:,5)==sid);
rid2 = find(r_all(:,:,6)==rid);
maxR(1,1) = lgrid(fid2);
maxR(1,2) = mgrid(fid2);
maxR(2,1) = lgrid(pid2);
maxR(2,2) = mgrid(pid2);
maxR(3,1) = lgrid(did2);
maxR(3,2) = mgrid(did2);
maxR(4,1) = lgrid(aid2);
maxR(4,2) = mgrid(aid2);
maxR(5,1) = lgrid(sid2);
maxR(5,2) = mgrid(sid2);
maxR(6,1) = lgrid(rid2);
maxR(6,2) = mgrid(rid2);

fie = squeeze(min(min(rmse_all(:,:,2),[],1),[],2));
pie = squeeze(min(min(rmse_all(:,:,3),[],1),[],2));
die = squeeze(min(min(rmse_all(:,:,4),[],1),[],2));
aie = squeeze(min(min(rmse_all(:,:,1),[],1),[],2));
sie = squeeze(min(min(rmse_all(:,:,5),[],1),[],2));
rie = squeeze(min(min(rmse_all(:,:,6),[],1),[],2));
fie2 = find(rmse_all(:,:,2)==fie);
pie2 = find(rmse_all(:,:,3)==pie);
die2 = find(rmse_all(:,:,4)==die);
aie2 = find(rmse_all(:,:,1)==aie);
sie2 = find(rmse_all(:,:,5)==sie);
rie2 = find(rmse_all(:,:,6)==rie);
minE(1,1) = lgrid(fie2);
minE(1,2) = mgrid(fie2);
minE(2,1) = lgrid(pie2);
minE(2,2) = mgrid(pie2);
minE(3,1) = lgrid(die2);
minE(3,2) = mgrid(die2);
minE(4,1) = lgrid(aie2);
minE(4,2) = mgrid(aie2);
minE(5,1) = lgrid(sie2);
minE(5,2) = mgrid(sie2);
minE(6,1) = lgrid(rie2);
minE(6,2) = mgrid(rie2);

