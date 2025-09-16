%% Area of locations

clear 
close all

%%
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t','area_t','ht');
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
area_t = double(area_t);
ht = double(ht);

%% fix lon shift
id=find(geolon_t<-180);
geolon_t(id)=geolon_t(id)+360;

%% 
surf_tmask = ht;
surf_tmask(~isnan(surf_tmask(:))) = ones;
surf_tmask(isnan(surf_tmask(:))) = zeros;

[nlon,nlat] = size(geolon_t);
t1 = [1:nlon*nlat]';
t1(:,2) = geolon_t(:);
t1(:,3) = geolat_t(:);
t1(:,4) = surf_tmask(:);
t1(:,5) = area_t(:);
t1(:,6) = ht(:);

ocean = find(t1(:,4)==1);
t2 = t1(ocean,:);

%% Georges Bank (Northeast Peak - 41 deg 43.92' N x 66 deg 32.18' W
%Southern Flank - 40 deg 57.95' N x 67 deg 18.91' W)
lon=find(t2(:,2)<=-66 & t2(:,2)>=-67);
lat=find(t2(:,3)<=42 & t2(:,3)>=41);
gid=intersect(lon,lat);

% W Scotian Shelf (42.4910,-65.4670)
lon=find(t2(:,2)<=-65 & t2(:,2)>=-66);
lat=find(t2(:,3)<=43 & t2(:,3)>=42);
wssid=intersect(lon,lat);

% C Scotian Shelf (43.5120,-62.4780)
lon=find(t2(:,2)<=-62 & t2(:,2)>=-63);
lat=find(t2(:,3)<=44 & t2(:,3)>=43);
cssid=intersect(lon,lat);

% E Scotian Shelf (45.0730,-59.1160)
lon=find(t2(:,2)<=-59 & t2(:,2)>=-60);
lat=find(t2(:,3)<=46 & t2(:,3)>=45);
essid=intersect(lon,lat);

% Greenland Sea (77.8710° N, 5.6501° W)
lon=find(t2(:,2)<=-5 & t2(:,2)>=-6);
lat=find(t2(:,3)<=78 & t2(:,3)>=77);
gsid=intersect(lon,lat);

% North Sea
lon=find(t2(:,2)<=4 & t2(:,2)>=3);
lat=find(t2(:,3)<=57 & t2(:,3)>=56);
nid=intersect(lon,lat);

% Norwegian Sea (68.8774° N, 3.1397° E)
lon=find(t2(:,2)<=4 & t2(:,2)>=3);
lat=find(t2(:,3)<=69 & t2(:,3)>=68);
nwid=intersect(lon,lat);

% Barents Sea (74.9884° N, 37.1064° E)
lon=find(t2(:,2)<=38 & t2(:,2)>=37);
lat=find(t2(:,3)<=75 & t2(:,3)>=74);
bsid=intersect(lon,lat);

% Eastern Bering Sea (M2 Southeastern Bering Sea (56.87°N, -164.06°W))
lon=find(t2(:,2)<=-164 & t2(:,2)>=-165);
lat=find(t2(:,3)<=57 & t2(:,3)>=56);
eid=intersect(lon,lat);

% Subarctic Pacific Gyre (Ocean Station Papa)
lon=find(t2(:,2)<=-145 & t2(:,2)>=-146);
lat=find(t2(:,3)<=51 & t2(:,3)>=50);
pid=intersect(lon,lat);

% HOT
lon=find(t2(:,2)<=-157 & t2(:,2)>=-158);
lat=find(t2(:,3)<=23 & t2(:,3)>=22);
hid=intersect(lon,lat);

% BATS
lon=find(t2(:,2)<=-64 & t2(:,2)>=-65);
lat=find(t2(:,3)<=32 & t2(:,3)>=31);
bid=intersect(lon,lat);

% Eastern Equatorial Pacific
lon=find(t2(:,2)<=-110 & t2(:,2)>=-111);
lat=find(t2(:,3)<=5.5 & t2(:,3)>=5);
qid=intersect(lon,lat);

% Peru Upwelling
lon=find(t2(:,2)<=-79 & t2(:,2)>=-80);
lat=find(t2(:,3)<=-12 & t2(:,3)>=-13);
uid=intersect(lon,lat);

%Subpolar W Pac station K2: 47oN, 160oE
lon=find(t2(:,2)<=161 & t2(:,2)>=160);
lat=find(t2(:,3)<=48 & t2(:,3)>=47);
kid=intersect(lon,lat);

%Subtropical W Pac station S1: 30oN, 145oE
lon=find(t2(:,2)<=146 & t2(:,2)>=145);
lat=find(t2(:,3)<=31 & t2(:,3)>=30);
sid=intersect(lon,lat);

%Southern Ocean: -60oN, 0.5oE
lon=find(t2(:,2)<=1 & t2(:,2)>=0);
lat=find(t2(:,3)<=-60 & t2(:,3)>=-61);
soid=intersect(lon,lat);

%S Atlantic: -20oN, -30oE
lon=find(t2(:,2)<=-30 & t2(:,2)>=-31);
lat=find(t2(:,3)<=-20 & t2(:,3)>=-20.5);
atlid=intersect(lon,lat);

%Indian ocean: -25oN, 105oE
lon=find(t2(:,2)<=106 & t2(:,2)>=105);
lat=find(t2(:,3)<=-25 & t2(:,3)>=-26);
indid=intersect(lon,lat);

%Pac Arctic: 80oN, -150oE
lon=find(t2(:,2)<=-150 & t2(:,2)>=-151);
lat=find(t2(:,3)<=81 & t2(:,3)>=80);
arcid=intersect(lon,lat);


%% Save
names={'Georges Bank','W Scotian Shelf','C Scotian Shelf','E Scotian Shelf',...
    'Greenland Sea','North Sea','Norwegian Sea','Barents Sea',...
    'Ocean Station Papa','Eastern Bering Sea','K2','S1',...
    'HOT','BATS','Eastern Equatorial Pacific','Peru Upwell',...
    'Southern Ocean','S Atl','Indian','Arctic'};
abbrev = {'GB','WSS','CSS','ESS',...
    'GS','NS','NwS','BS',...
    'OSP','EBS','K2','S1',...
    'HOT','BATS','EEP','PUp',...
    'SO','SAtl','Ind','Arc'};

ids(1,1)=gid;
ids(2,1)=wssid;
ids(3,1)=cssid;
ids(4,1)=essid;
ids(5,1)=gsid;
ids(6,1)=nid;
ids(7,1)=nwid;
ids(8,1)=bsid;
ids(9,1)=pid;
ids(10,1)=eid;
ids(11,1)=kid;
ids(12,1)=sid;
ids(13,1)=hid;
ids(14,1)=bid;
ids(15,1)=qid;
ids(16,1)=uid;
ids(17,1)=soid;
ids(18,1)=atlid;
ids(19,1)=indid;
ids(20,1)=arcid;

lons = t2(ids,2);
lats = t2(ids,3);
area = t2(ids,5);
depth = t2(ids,6);

%%
T=table(names',ids,lons,lats,area,depth,...
    'VariableNames',{'Location','ID','Lon','Lat','Area','Depth'});
writetable(T,'/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/core_grid_360x200_id_locs_area_dep.csv','Delimiter',',');
save('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/core_grid_360x200_id_locs_area_dep.mat','T','depth','ids','abbrev');

%%
glocs = zeros(ni,nj);
glocs(GRD.ID(ids)) = ones;

%%
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,glocs)
for n=1:length(ids)
    textm(geolat_t(GRD.ID(ids(n))),geolon_t(GRD.ID(ids(n))),abbrev{n});
    hold on;
end
cmocean('dense')
load coastlines;                     
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1]);
hcb = colorbar('h');
print('-dpng',[ppath 'map_CORE_locs.png'])
