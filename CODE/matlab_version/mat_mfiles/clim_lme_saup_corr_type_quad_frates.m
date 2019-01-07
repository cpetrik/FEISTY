%POEM catch vs. SAUP catch by LME
%Find fishing selectivity that gives highest corr for each type

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(3,:);
cmap4(2,:)=cmap_ppt(1,:);
cmap4(3,:)=cmap_ppt(5,:);

AREA_OCN = max(area,1);

lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = [dp lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%%
Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1956-2005 SAUP average
id = find(yr>1955 & yr<=2005);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));
%Top 10 yrs SAUP
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
slme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);

keep = [1:9,11:34,36:38,46:54,59:60,62,65];

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% Loop over POEM params
fqs = [0.25 0.5:0.5:3.5];
pqs = 0.25:0.25:1.5;
dqs = 1:7;
frate = 0.1;

rall = NaN*ones(length(fqs),length(pqs),length(dqs));
rF = rall;
rP = rall;
rD = rall;
rmse = rall;
rmseF = rall;
rmseP = rall;
rmseD = rall;

sname = 'Clim_';
dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
dpath = [datap char(dp) '/'];
ppath = [pp char(dp) '/'];

for fq=1:length(fqs)
    for pq=1:length(pqs)
        for dq=1:length(dqs)
            MFsel = fqs(fq);
            LPsel = pqs(pq);
            LDsel = dqs(dq);
            tF = num2str(1000+int64(100*frate*MFsel));
            tP = num2str(1000+int64(100*frate*LPsel));
            tD = num2str(1000+int64(100*frate*LDsel));
            
            charv = ['fish_F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
            harv = ['F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
            cfile = [harv,'_cm20_m-b175-k09_D075_J100_A050_Sm025_nmort1_BE05_CC100_RE00100'];

            if exist([dpath 'LME_clim_fished_',harv,'_loop_' dp '.mat'],'file')
            load([dpath 'LME_clim_fished_',harv,'_loop_' dp '.mat'],'lme_mcatch');
            
            close all 
            
            %% POEM LME biomass in MT
            plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
            plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
            plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
            plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
            % MT/km2
            plme_mcatch = plme_mcatch ./ lme_area_km2;
            plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
            plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
            plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
            
            %log10 Difference
            diff_catch = (log10(plme_mcatch) - log10(slme_mcatch10)) ./ log10((slme_mcatch10+eps));
            Fdiff_catch = (log10(plme_Fmcatch) - log10(Flme_mcatch10)) ./ log10((Flme_mcatch10+eps));
            Pdiff_catch = (log10(plme_Pmcatch) - log10(Plme_mcatch10)) ./ log10((Plme_mcatch10+eps));
            Ddiff_catch = (log10(plme_Dmcatch) - log10(Dlme_mcatch10)) ./ log10((Dlme_mcatch10+eps));
            
            code = [1:66]';
            T = table(code,slme_mcatch10,plme_mcatch,Flme_mcatch10,plme_Fmcatch,...
                Plme_mcatch10,plme_Pmcatch,Dlme_mcatch10,plme_Dmcatch,'VariableNames',...
                {'lme','saupAll','poemAll','saupF','poemF','saupP','poemP','saupD','poemD'});
            writetable(T,[dpath 'LME_saup_catch_clim_fished',harv,'.csv']);
            
            %% Plots
            l10s=log10(slme_mcatch10+eps);
            l10p=log10(plme_mcatch);
            l10sF=log10(Flme_mcatch10+eps);
            l10pF=log10(plme_Fmcatch);
            l10sP=log10(Plme_mcatch10+eps);
            l10pP=log10(plme_Pmcatch);
            l10sD=log10(Dlme_mcatch10+eps);
            l10pD=log10(plme_Dmcatch);
            
            
            %% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
             
            % Stats
            %r
            rall(fq,pq,dq)=corr(l10s(keep),l10p(keep));
            rF(fq,pq,dq)=corr(l10sF(keep),l10pF(keep));
            rP(fq,pq,dq)=corr(l10sP(keep),l10pP(keep));
            rD(fq,pq,dq)=corr(l10sD(keep),l10pD(keep));
            
            %root mean square error
            o=l10s(keep);
            p=l10p(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmse(fq,pq,dq) = sqrt(num/n);
            
            o=l10sF(keep);
            p=l10pF(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmseF(fq,pq,dq) = sqrt(num/n);
            
            o=l10sP(keep);
            p=l10pP(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmseP(fq,pq,dq) = sqrt(num/n);
            
            o=l10sD(keep);
            p=l10pD(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmseD(fq,pq,dq) = sqrt(num/n);
            
            end %if
          
        end
    end
end

%% Plots of r and RMSE
cfile2 = ['_cm20_m-b175-k09_D075_J100_A050_Sm025_nmort1_BE05_CC100_RE00100'];

fqs2 = [0.25 0.5:0.5:4];
pqs2 = 0.25:0.25:1.75;
dqs2 = 1:8;
[pgrid1,fgrid1]=meshgrid(pqs2,fqs2);
[dgrid1,fgrid2]=meshgrid(dqs2,fqs2);
[dgrid2,pgrid2]=meshgrid(dqs2,pqs2);
nf = length(fqs);
np = length(pqs);
nd = length(dqs);

Rall= NaN*ones(nf+1,np+1,nd+1);
RF  = NaN*ones(nf+1,np+1,nd+1);
RP  = NaN*ones(nf+1,np+1,nd+1);
RD  = NaN*ones(nf+1,np+1,nd+1);
Rall(1:nf,1:np,1:nd)=rall;
RF(1:nf,1:np,1:nd)=rF;
RP(1:nf,1:np,1:nd)=rP;
RD(1:nf,1:np,1:nd)=rD;

% for j=1:nd
%     figure(2)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,Rall(:,:,j))
%     colorbar
%     caxis([0.3 0.5])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('Corr all fish')
%     %%print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrAll_frates_FP'])
%     
%     figure(3)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,RF(:,:,j))
%     colorbar
%     caxis([0.3 0.45])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('Corr F')
%     %%print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrF_frates_FP'])
%     
%     figure(4)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,RP(:,:,j))
%     colorbar
%     caxis([0.5 0.6])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('Corr P')
%     %%print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrP_frates_FP'])
%     
%     figure(5)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,RD(:,:,j))
%     colorbar
%     caxis([0.5 0.55])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('Corr D')
%     stamp(cfile2)
%     %%print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_FP'])
% end
%%
for j=1:np
    figure(6)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(Rall(:,j,:)))
    colorbar
    caxis([0.3 0.5])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr all fish')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrAll_frates_FD'])
    
    figure(7)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RF(:,j,:)))
    colorbar
    caxis([0.3 0.45])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr F')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrF_frates_FD'])
    
    figure(8)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RP(:,j,:)))
    colorbar
    caxis([0.5 0.6])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr P')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrP_frates_FD'])
    
    figure(9)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RD(:,j,:)))
    colorbar
    caxis([0.5 0.55])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr D')
    stamp(cfile2)
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_FD'])
end
%%
for j=1:nf
    figure(10)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(Rall(j,:,:)))
    colorbar
    caxis([0.3 0.5])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr all fish')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrAll_frates_PD'])
    
    figure(11)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RF(j,:,:)))
    colorbar
    caxis([0.3 0.45])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr F')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrF_frates_PD'])
    
    figure(12)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RP(j,:,:)))
    colorbar
    caxis([0.5 0.6])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr P')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrP_frates_PD'])
    
    figure(13)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RD(j,:,:)))
    colorbar
    caxis([0.5 0.55])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr D')
    stamp(cfile2)
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_PD'])
    
end

%%
Rmseall  = NaN*ones(nf+1,np+1,nd+1);
RmseF  = NaN*ones(nf+1,np+1,nd+1);
RmseP  = NaN*ones(nf+1,np+1,nd+1);
RmseD  = NaN*ones(nf+1,np+1,nd+1);
Rmseall(1:nf,1:np,1:nd)=rmse;
RmseF(1:nf,1:np,1:nd)=rmseF;
RmseP(1:nf,1:np,1:nd)=rmseP;
RmseD(1:nf,1:np,1:nd)=rmseD;

% for j=1:nd
%     figure(14)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,Rmseall(:,:,j))
%     colorbar
%     caxis([0.525 0.675])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('RMSE all fish')
%     %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseAll_frates_FP'])
%     
%     figure(15)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,RmseF(:,:,j))
%     colorbar
%     caxis([1.15 1.25])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('RMSE F')
%     %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseF_frates_FP'])
%     
%     figure(16)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,RmseP(:,:,j))
%     colorbar
%     caxis([1 1.1])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('RMSE P')
%     %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseP_frates_FP'])
%     
%     figure(17)
%     subplot(2,3,j)
%     pcolor(fgrid1,pgrid1,RmseD(:,:,j))
%     colorbar
%     caxis([0.6 1])
%     %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
%     xlabel('F sel')
%     ylabel('P sel')
%     title('RMSE D')
%     stamp(cfile2)
%     %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_FP'])
% end
for j=1:np
    figure(18)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(Rmseall(:,j,:)))
    colorbar
    caxis([0.525 0.675])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE all fish')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseAll_frates_FD'])
    
    figure(19)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RmseF(:,j,:)))
    colorbar
    caxis([1.15 1.25])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE F')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseF_frates_FD'])
    
    figure(20)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RmseP(:,j,:)))
    colorbar
    caxis([1 1.1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE P')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseP_frates_FD'])
    
    figure(21)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RmseD(:,j,:)))
    colorbar
    caxis([0.6 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE D')
    stamp(cfile2)
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_FD'])
end
%%
for j=1:nf
    figure(22)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(Rmseall(j,:,:)))
    colorbar
    caxis([0.525 0.675])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE all fish')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseAll_frates_PD'])
    
    figure(23)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RmseF(j,:,:)))
    colorbar
    caxis([1.15 1.25])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE F')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseF_frates_PD'])
    
    figure(24)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RmseP(j,:,:)))
    colorbar
    caxis([1 1.1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE P')
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseP_frates_PD'])
    
    figure(25)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RmseD(j,:,:)))
    colorbar
    caxis([0.6 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE D')
    stamp(cfile2)
    %print('-dpng',[ppath 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_PD'])
    
end

%%
[pgrid3,fgrid3,dgrid3]=meshgrid(pqs2,fqs2,dqs2);

vals = NaN*ones(2,6);
qr = NaN*ones(3,6);
qrmse = NaN*ones(3,6);

corr_all = 0.25*(normalize(Rall)+normalize(RF)+normalize(RP)+normalize(RD));
rmse_all = normalize(Rmseall)+normalize(RmseF)+normalize(RmseP)+normalize(RmseD);
corr_PD = 0.5*(normalize(RP)+normalize(RD));
rmse_PD = normalize(RmseP)+normalize(RmseD);

ax1 = find(Rall(:)==max(Rall(:)));
fx1 = find(RF(:)==max(RF(:)));
px1 = find(RP(:)==max(RP(:)));
dx1 = find(RD(:)==max(RD(:)));
lx1 = find(corr_all(:)==max(corr_all(:)));
pdx1 = find(corr_PD(:)==max(corr_PD(:)));

ax2 = find(Rmseall(:)==min(Rmseall(:)));
fx2 = find(RmseF(:)==min(RmseF(:)));
px2 = find(RmseP(:)==min(RmseP(:)));
dx2 = find(RmseD(:)==min(RmseD(:)));
lx2 = find(rmse_all(:)==min(rmse_all(:)));
pdx2 = find(rmse_PD(:)==min(rmse_PD(:)));

vals(1,1) = Rall(ax1);
vals(1,2) = RF(fx1);
vals(1,3) = RP(px1);
vals(1,4) = RD(dx1);
vals(1,5) = corr_all(lx1);
vals(1,6) = corr_PD(pdx1);
vals(2,1) = Rmseall(ax2);
vals(2,2) = RmseF(fx2);
vals(2,3) = RmseP(px2);
vals(2,4) = RmseD(dx2);
vals(2,5) = rmse_all(lx2);
vals(2,6) = rmse_PD(pdx2);

qr(1,1) = fgrid3(ax1);
qr(1,2) = fgrid3(fx1);
qr(1,3) = fgrid3(px1);
qr(1,4) = fgrid3(dx1);
qr(1,5) = fgrid3(lx1);
qr(1,6) = fgrid3(pdx1);
qr(2,1) = pgrid3(ax1);
qr(2,2) = pgrid3(fx1);
qr(2,3) = pgrid3(px1);
qr(2,4) = pgrid3(dx1);
qr(2,5) = pgrid3(lx1);
qr(2,6) = pgrid3(pdx1);
qr(3,1) = dgrid3(ax1);
qr(3,2) = dgrid3(fx1);
qr(3,3) = dgrid3(px1);
qr(3,4) = dgrid3(dx1);
qr(3,5) = dgrid3(lx1);
qr(3,6) = dgrid3(pdx1);

qrmse(1,1) = fgrid3(ax2);
qrmse(1,2) = fgrid3(fx2);
qrmse(1,3) = fgrid3(px2);
qrmse(1,4) = fgrid3(dx2);
qrmse(1,5) = fgrid3(lx2);
qrmse(1,6) = fgrid3(pdx2);
qrmse(2,1) = pgrid3(ax2);
qrmse(2,2) = pgrid3(fx2);
qrmse(2,3) = pgrid3(px2);
qrmse(2,4) = pgrid3(dx2);
qrmse(2,5) = pgrid3(lx2);
qrmse(2,6) = pgrid3(pdx2);
qrmse(3,1) = dgrid3(ax2);
qrmse(3,2) = dgrid3(fx2);
qrmse(3,3) = dgrid3(px2);
qrmse(3,4) = dgrid3(dx2);
qrmse(3,5) = dgrid3(lx2);
qrmse(3,6) = dgrid3(pdx2);

qr(:,7) = mean(qr,2);
qrmse(:,7) = mean(qrmse,2);

%%
save([dpath 'Clim_fished_loop_LME_SAUP_catch_comp_qs.mat'],'rall','rF','rP',...
    'rD','rmse','rmseF','rmseP','rmseD','fgrid3','pgrid3','dgrid3',...
    'vals','qr','qrmse')
csvwrite([dpath 'Clim_fished_loop_LME_SAUP_catch_comp_qs_corr.csv'],qr);
csvwrite([dpath 'Clim_fished_loop_LME_SAUP_catch_comp_qs_rmse.csv'],qrmse);

