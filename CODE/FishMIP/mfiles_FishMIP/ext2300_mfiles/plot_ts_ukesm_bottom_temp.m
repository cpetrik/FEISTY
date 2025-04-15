% SSP 534-over 2015-2100
% plot theta-bot because looks off from tob in other scenarios

clear 
close all

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/WGs/2300/testing_forcing/';

%%
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp534over/';

load([fpath 'ukesm_ssp534_temp_btm_monthly_2040_2100.mat']);

temp_btm(temp_btm > 1.0e19) = nan;

ssp534_Tb = temp_btm;
ssp534_yr = yr;

clear temp_btm yr

%% 
spath='/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/ssp585/';

load([spath 'ukesm_ssp585_temp_btm_monthly_2015_2300.mat']);

ssp585_Tb = temp_btm;
ssp585_yr = yr;

clear temp_btm yr

load('/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/gridspec_ukesm_cmip6_2300.mat','deptho')

%% find deep & shallow locations to plot

wid = find(deptho >0);

[ni,nj] = size(deptho);
locs34 = nan*ones(ni,nj);
locs34(wid) = wid;

figure
pcolor(locs34); shading flat;
colorbar

[LAT,LON] = meshgrid(lat,lon);

figure
pcolor(LAT); shading flat;
colorbar
title('lat')

figure
pcolor(LON); shading flat;
colorbar
title('lon')

%% find some locs
% Georges Bank (Northeast Peak - 41 deg 43.92' N x 66 deg 32.18' W
%Southern Flank - 40 deg 57.95' N x 67 deg 18.91' W)
lon=find(LON(:)<=-67 & LON(:)>=-68);
lat=find(LAT(:)<=41.5 & LAT(:)>=41);
gid=intersect(lon,lat);

% C Scotian Shelf (43.5120,-62.4780)
lon=find(LON(:)<=-62 & LON(:)>=-63);
lat=find(LAT(:)<=44 & LAT(:)>=43.5);
cssid=intersect(lon,lat);

% Greenland Sea (77.8710° N, 5.6501° W)
lon=find(LON(:)<=-5.5 & LON(:)>=-6);
lat=find(LAT(:)<=78 & LAT(:)>=77.5);
gsid=intersect(lon,lat);

% North Sea
lon=find(LON(:)<=4 & LON(:)>=3.25);
lat=find(LAT(:)<=57 & LAT(:)>=56.5);
nid=intersect(lon,lat);

% Norwegian Sea (68.8774° N, 3.1397° E)
lon=find(LON(:)<=4 & LON(:)>=3);
lat=find(LAT(:)<=69 & LAT(:)>=68);
nwid=intersect(lon,lat);

% Barents Sea (74.9884° N, 37.1064° E)
lon=find(LON(:)<=38 & LON(:)>=37);
lat=find(LAT(:)<=75.4 & LAT(:)>=74.5);
bsid=intersect(lon,lat);

% Eastern Bering Sea (M2 Southeastern Bering Sea (56.87°N, -164.06°W))
lon=find(LON(:)<=-164 & LON(:)>=-165);
lat=find(LAT(:)<=57 & LAT(:)>=56.5);
eid=intersect(lon,lat);

% Subarctic Pacific Gyre (Ocean Station Papa)
lon=find(LON(:)<=-145 & LON(:)>=-146);
lat=find(LAT(:)<=51 & LAT(:)>=50);
pid=intersect(lon,lat);

% HOT (22° 45'N, 158° 00'W)
lon=find(LON(:)<=-157.35 & LON(:)>=-158);
lat=find(LAT(:)<=23 & LAT(:)>=22.45);
hid=intersect(lon,lat);

% BATS (31 50'N, 64 10'W)
lon=find(LON(:)<=-63.95 & LON(:)>=-65);
lat=find(LAT(:)<=31.9 & LAT(:)>=31.1);
bid=intersect(lon,lat);

% Eastern Equatorial Pacific
lon=find(LON(:)<=-110.3 & LON(:)>=-111);
lat=find(LAT(:)<=5.4 & LAT(:)>=4.5);
qid=intersect(lon,lat);

% Peru Upwelling
lon=find(LON(:)<=-79 & LON(:)>=-80);
lat=find(LAT(:)<=-12.3 & LAT(:)>=-13);
uid=intersect(lon,lat);

%Subpolar W Pac station K2: 47oN, 160oE
lon=find(LON(:)<=161 & LON(:)>=160);
lat=find(LAT(:)<=47.5 & LAT(:)>=47);
kid=intersect(lon,lat);

%Subtropical W Pac station S1: 30oN, 145oE
lon=find(LON(:)<=146 & LON(:)>=145);
lat=find(LAT(:)<=30.5 & LAT(:)>=30);
sid=intersect(lon,lat);


%% Save
names={'Georges Bank','C Scotian Shelf',...
    'Greenland Sea','North Sea','Norwegian Sea','Barents Sea',...
    'Ocean Station Papa','Eastern Bering Sea','K2','S1',...
    'HOT','BATS','Eastern Equatorial Pacific','Peru Upwell'};
abbrev = {'GB','CSS',...
    'GS','NS','NwS','BS',...
    'OSP','EBS','K2','S1',...
    'HOT','BATS','EEP','PUp'};

ids(1,1)=gid;
ids(2,1)=cssid;
ids(3,1)=gsid;
ids(4,1)=nid;
ids(5,1)=nwid;
ids(6,1)=bsid;
ids(7,1)=pid;
ids(8,1)=eid;
ids(9,1)=kid;
ids(10,1)=sid;
ids(11,1)=hid;
ids(12,1)=bid;
ids(13,1)=qid;
ids(14,1)=uid;

lons = LON(ids);
lats = LAT(ids);
depth = deptho(ids);

ids(:,2) = lons;
ids(:,3) = lats;
ids(:,4) = depth;

T=array2table(ids,'RowNames',names,...
    'VariableNames',{'ID','Lon','Lat','Depth'});

wpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/CMIP6/UKESM1-0-LL/';
writetable(T,[wpath 'cesm_grid_id_locs_area_dep.csv'],'Delimiter',',');
save([wpath 'cesm_grid_id_locs_dep.mat'],'T','ids','abbrev');

%% Temps at specific locs
[ni,nj,nt] = size(ssp534_Tb);
[nu,nv,nw] = size(ssp585_Tb);

ssp534_Tb = reshape(ssp534_Tb, ni*nj,nt);
ssp585_Tb = reshape(ssp585_Tb, nu*nv,nw);

%%
figure
subplot(2,2,1)
plot(ssp585_yr(1:12:end),ssp585_Tb(nid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(nid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - North Sea ' num2str(round(deptho(nid))) 'm'])
legend('585','534')
legend('location','northwest')

subplot(2,2,2)
plot(ssp585_yr(1:12:end),ssp585_Tb(eid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(eid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - E Bering Sea ' num2str(round(deptho(eid))) 'm'])

subplot(2,2,3)
plot(ssp585_yr(1:12:end),ssp585_Tb(kid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(kid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - K2 ' num2str(round(deptho(kid))) 'm'])

subplot(2,2,4)
plot(ssp585_yr(1:12:end),ssp585_Tb(bid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(bid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - BATS ' num2str(round(deptho(bid))) 'm'])

print('-dpng',[pp 'UKESM_bottom_temp_locs.png'])


%%
figure
subplot(2,2,1)
plot(ssp585_yr(1:12:end),ssp585_Tb(gid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(gid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - Greenland Sea ' num2str(round(deptho(gid))) 'm'])
legend('585','534')
legend('location','northwest')

subplot(2,2,2)
plot(ssp585_yr(1:12:end),ssp585_Tb(uid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(uid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - Peru Upwelling ' num2str(round(deptho(uid))) 'm'])

subplot(2,2,3)
plot(ssp585_yr(1:12:end),ssp585_Tb(sid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(sid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - S1 ' num2str(round(deptho(sid))) 'm'])

subplot(2,2,4)
plot(ssp585_yr(1:12:end),ssp585_Tb(hid,1:12:end),'r','LineWidth',2); hold on
plot(ssp534_yr(1:12:end),ssp534_Tb(hid,1:12:end),'color',[0 0.75 0.5],'LineWidth',2);
title(['UKESM Tbtm - HOT ' num2str(round(deptho(hid))) 'm'])

print('-dpng',[pp 'UKESM_bottom_temp_locs2.png'])


