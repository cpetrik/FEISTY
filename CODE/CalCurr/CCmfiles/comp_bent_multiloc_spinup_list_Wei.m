% Calculate different skill metrics for global bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

% np1 = 'Dc_TrefO_Hold_cmax-metab_MFeqMP_fcrit40_D100_nmort2_BE05_CC050_RE1000';
% np2 = 'Dc_TrefO_Hold_cmax-metab_MFeqMP_fcrit40_D100_nmort2_BE05_CC050_RE0500';
% np3 = 'Dc_TrefO_Hold_cmax-metab_MFeqMP_fcrit40_D100_nmort2_BE05_CC050_RE0100';
% np4 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0500';
% np5 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0100';
% np6 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0050';
% np7 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0010';
% np8 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE10_CC050_RE0100';
% np9 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE10_CC050_RE0050';
% np10 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC100_RE0100';
% dp = {np1;np2;np3;np4;np5;np6;np7;np8;np9;np10};
%
% sims = {'212-05-05-10';'212-05-05-05';'212-05-05-01';'443-05-05-05';'443-05-05-01';'443-05-05-005';...
%     '443-05-05-001';'443-10-05-01';'443-10-05-005';'443-05-10-01'};

% np1 = 'Dc_enc70_cmax-metab20_b175_k08_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
% np2 = 'Dc_enc70_cmax-metab20_b180_k08_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
% np3 = 'Dc_enc70_cmax-metab20_b175_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
% np4 = 'Dc_enc70_cmax-metab20_b180_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
np1 = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
np2 = 'Dc_enc70_cmax-metab20_b175_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
np3 = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE10_CC050_lgRE00100_mdRE00100';
np4 = 'Dc_enc70_cmax-metab20_b175_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE10_CC100_lgRE00100_mdRE00100';
dp = {np1;np2;np3;np4};

% sims = {'175-08';'180-08';'175-09';'180-09'};
sims = {'05-05';'05-10';'10-05';'10-10'};

%% Wei benthic biomass
seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei2010_Global_seafloor_biomass.csv',1,0);
Wcol = {'latitude','longitude','depth','bact.biom.mean','meio.biom.mean',...
    'macro.biom.mean','mega.biom.mean','inv.biom.mean','fis.biom.mean'};
Wcol = Wcol';

% all mean biomasses in log10 mg C/m2
invert = seafl(:,8);
fish = seafl(:,9);
% convert to g WW/m2
invert = 10.^(invert) * 1e-3 * 9.0;
fish = 10.^(fish) * 1e-3 * 9.0;

%% put on same grid as POEM output
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

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
test=seafl(:,2);
id=find(test>80);
test(id)=test(id)-360;

geolat_c = double(lat);
geolon_c = double(lon) - 280;

Zi = griddata(seafl(:,1),test,invert,geolat_c,geolon_c);
Zf = griddata(seafl(:,1),test,fish,geolat_c,geolon_c);
% Zi = Zi.*lmask;
% Zf = Zf.*lmask;

%% ocean cells
% ocean=zeros(ni,nj);
% ocean(grid(:,1))=ones(size(grid(:,1)));
% figure(55)
% surf(geolon_t,geolat_t,ocean); view(2); hold on;

%% Maps
% Inv
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_c,geolon_c,log10(Zi))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.5 0.5]);
hcb = colorbar('h');
ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Wei log10 mean benthic invert biomass (g m^-^2)')
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_inverts_Robinson.png')


% Fish
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_c,geolon_c,log10(Zf))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Wei log10 mean benthic fish biomass (g m^-^2)')
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_fish_Robinson.png')


%% Skill metrics

% Stowe et al 2009
% Univariate (model vs observ)
% 1. correlation coefficient [-1 1]; want close to 1
% 2. root mean square error; want close to 0; considering the magnitude rather than the direction
% 3. average error (bias); want close to 0; neg & pos discrepancies can cancel each other
% 4. average absolute error; want close to 0; considering the magnitude rather than the direction
% 5. modeling efficiency; want close to 1; <0 means avg obs better than model; how well a model predicts relative to avg obs

metrics={'r','RMSE','AE','AAE','MEF'};
metrics=metrics';

%Sometimes it is appropriate to log-transform the observations and predictions before
%calculating goodness-of-fit statistics so that differences between predicted and observed
%values will not be highly skewed and dominated by a small proportion of high values.

r = NaN*ones(length(dp),1);
RMSE = r;
AE = r;
AAE = r;
MEF = r;
NSTD = r;
RMSD = r;

%%
iskill=NaN*ones(7,length(dp));
fskill=NaN*ones(7,length(dp));
ibest=cell(5,1);
fbest=cell(5,1);

for i=1:length(dp)
    close all
    
    %%
    cfile = dp{i};
    fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
    ppath = [pp cfile];
    load([fpath 'Means_Climatol_All_fish03_' cfile '.mat'],'md_mean','ld_mean','b_mean');
    
    %% Model bent & dem
    Zmd=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);
    
    Zmd(ID)=md_mean;
    Zld(ID)=ld_mean;
    Zb(ID)=b_mean;
    
    Zd = Zmd+Zld;
    
    % Bent
%     figure(3)
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%     surfm(geolat_t,geolon_t,log10(Zb))
%     colormap('jet')
%     load coast;                     %decent looking coastlines
%     h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     caxis([-2.5 0.5]);
%     hcb = colorbar('h');
%     ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
%     set(gcf,'renderer','painters')
%     title('Climatology log10 mean Benthic inverts (g m^-^2)')
%     stamp(cfile)
%     print('-dpng',[ppath 'Climatol_All_fish03_global_BENT.png'])
    
    % D
    figure(4)
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(geolat_t,geolon_t,log10(Zd))
    colormap('jet')
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    caxis([-1 1]);
    hcb = colorbar('h');
    ylim(hcb,[-1 1])                   %Set color axis if needed
    set(gcf,'renderer','painters')
    title('Climatology log10 mean M&L D (g m^-^2)')
    stamp(cfile)
    print('-dpng',[ppath 'Climatol_All_fish03_global_Demersal.png'])
    
    
    %% Inverts
    obs = real(log(Zi(:)));
    model = real(log(Zb(:)));
    
    n = length(obs);
    
    o=obs;
    p=model;
    %             o(isnan(o))=0;
    %             p(isnan(p))=0;
    
    omean=repmat(nanmean(o),n,1);
    pmean=repmat(nanmean(p),n,1);
    osig=nanstd(o);
    psig=nanstd(p);
    
    % corr coeff
    num=nansum((o-omean).*(p-pmean));
    d1=nansum((o-omean).^2);
    d2=nansum((p-pmean).^2);
    den=sqrt(d1*d2);
    iskill(1,i) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    iskill(2,i) = sqrt(num/n);
    
    % average error
    iskill(3,i) = nansum(p-o) / n;
    
    % average absolute error
    iskill(4,i) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    iskill(5,i) = (num1-num2)/num1;
    
    % Taylor normalized std
    iskill(6,i) = psig/osig;
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    iskill(7,i) = sqrt(q1/n);
    
    
    %% Fish
    obs = real(log(Zf(:)));
    model = real(log(Zd(:)));
    
    n = length(obs);
    
    o=obs;
    p=model;
    
    omean=repmat(nanmean(o),n,1);
    pmean=repmat(nanmean(p),n,1);
    osig=nanstd(o);
    psig=nanstd(p);
    
    % corr coeff
    num=nansum((o-omean).*(p-pmean));
    d1=nansum((o-omean).^2);
    d2=nansum((p-pmean).^2);
    den=sqrt(d1*d2);
    fskill(1,i) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    fskill(2,i) = sqrt(num/n);
    
    % average error
    fskill(3,i) = nansum(p-o) / n;
    
    % average absolute error
    fskill(4,i) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    fskill(5,i) = (num1-num2)/num1;
    
    % Taylor normalized std
    fskill(6,i) = psig/osig;
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    fskill(7,i) = sqrt(q1/n);
    
end %sims

%% Plot results
cfile2 = 'Dc_enc70_cmax-metab20_b175_k09_fcrit20_D075_J100_A050_Sm025_nmort1_lgRE00100_mdRE00100_BECCtests';

% Bar graphs
figure(5)
subplot(2,2,1)
bar(iskill(1,:),'k')
ylabel('Correlation coefficient')
ylim([0 0.4])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.01,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
stamp(cfile2)
title('Inverts skill','HorizontalAlignment','left')

subplot(2,2,2)
bar(iskill(2,:),'k')
ylabel('Root mean square error')
ylim([0 2])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,3)
bar(iskill(3,:),'k')
ylabel('Average error')
ylim([0 1])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,4)
bar(iskill(5,:),'k')
ylabel('Modeling efficiency')
ylim([-0.8 0])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.82,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
print('-dpng',[pp cfile2 '_Clim_skill_Wei_inverts'])

%%
figure(6)
subplot(2,2,1)
bar(fskill(1,:),'k')
ylabel('Correlation coefficient')
ylim([0 0.05])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.001,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
stamp(cfile2)
title('Fish skill','HorizontalAlignment','left')

subplot(2,2,2)
bar(fskill(2,:),'k')
ylabel('Root mean square error')
ylim([0 2.5])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,3)
bar(fskill(3,:),'k')
ylabel('Average error')
ylim([0 1.5])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,4)
bar(fskill(5,:),'k')
ylabel('Modeling efficiency')
ylim([-1.5 0])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-1.55,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
print('-dpng',[pp cfile2 '_Clim_skill_Wei_fish'])

%% Taylor diagrams
% % INVERTS
% [rmsd,it]=sort(iskill(7,:),'descend');
% theta=acos(iskill(1,it));    %corr coeff
% rho=iskill(6,it);            %stdev
% simtex=sims(it);
% simtex{length(sims)+1}='obs';
% 
% % % Just those <2
% % id = find(rho<2);
% % theta = theta(id);
% % rho = rho(id);
% % simtex=sims(id);
% % simtex{length(id)+1}='obs';
% 
% tr=0;
% rr=1;
% figure(7)
% h0=polar(0,1.5,'.'); hold on;
% set(h0,'color','w');
% for s=1:length(dp)
%     h=polar(theta(s),rho(s),'.'); hold on;
%     set(h,'color',cm{s},'MarkerSize',25);
%     %     set(h,'MarkerSize',25);
% end
% h2=polar(tr,rr,'k*');
% set(h2,'MarkerSize',10);
% axis([0 1.5 0 1.5])
% title('Inverts Taylor diagram')
% legend([' '; simtex])
% legend('location','northeast')
% print('-dpng',[pp cfile2 '_Clim_global_invert_Taylor_Wei'])
% 
% %% FISH
% 
% [rmsd,it]=sort(fskill(7,:),'descend');
% theta=acos(fskill(1,it));    %corr coeff
% rho=fskill(6,it);            %stdev
% simtex=sims(it);
% simtex{length(sims)+1}='obs';
% 
% % Just those <1.5
% % id = find(rho<1.5);
% % theta = theta(id);
% % rho = rho(id);
% % simtex=sims(id);
% % simtex{length(id)+1}='obs';
% 
% tr=0;
% rr=1;
% figure(8)
% h0=polar(0,1.5,'.'); hold on;
% set(h0,'color','w');
% for s=1:length(dp)
%     h=polar(theta(s),rho(s),'.'); hold on;
%     set(h,'color',cm{s},'MarkerSize',25);
%     %     set(h,'MarkerSize',25);
% end
% h2=polar(tr,rr,'k*');
% set(h2,'MarkerSize',10);
% axis([0 1.5 0 1.5])
% title('Fish Taylor diagram')
% legend([' '; simtex])
% legend('location','northeast')
% print('-dpng',[pp cfile2 '_Clim_global_fish_Taylor_Wei'])

%% Best
one=[1;5];
zer=[2;3;4];
for z=1:5
    Z1=sum(z==one);
    Z2=sum(z==zer);
    if (Z1>0)
        ibest{z,1}=sims(find(min(abs(1-(iskill(z,:))))==abs(1-(iskill(z,:)))));
        fbest{z,1}=sims(find(min(abs(1-(fskill(z,:))))==abs(1-(fskill(z,:)))));
    elseif (Z2>0)
        ibest{z,1}=sims(find(min(abs(0-(iskill(z,:))))==abs(0-(iskill(z,:)))));
        fbest{z,1}=sims(find(min(abs(0-(fskill(z,:))))==abs(0-(fskill(z,:)))));
    end
end

Ti=table(metrics,ibest);
Tf=table(metrics,fbest);

%%
save(['/Volumes/GFDL/NC/Matlab_new_size/' cfile2 '_global_skill_Wei.mat'],...
    'invert','fish','metrics','iskill','fskill','Ti','Tf');



