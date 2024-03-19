% Put fishing mortality onto grid

clear
close all

%% grid info
Cdir = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';

% Depth, lat, lon
load([Cdir 'feisty_ipsl_gridspec.mat']);
load([Cdir 'Data_grid_nemuro_ipsl.mat']);

[ni,nj] = size(LON);

WID = GRD.ID;
NID = GRD.N;

clatlim=[25 50];
clonlim=[-140 -110];
load coastlines;

%% temp scaling
load([Cdir 'feisty_ipsl_tp_1980-2100.mat']);
load([Cdir 'feisty_ipsl_tb_1980-2100.mat']);

%% take mean temp over hist time period = 1980-2010
% time
year = 1980:(1/12):2010;
yrall = 1980:(1/12):2101;
[y,yid] = intersect(year,yrall);

mtp = nanmean(TEMP_AVG_200M(:,:,yid),3);
mtb = nanmean(TEMP_BOT(:,:,yid),3);

vmtp = mtp(WID);
vmtb = mtb(WID);

%% load f/fmsy
alt = {'assessment'};

i=1;
alt1 = alt{i};
spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
fpath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
    'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
load([fpath 'grid_mortality_all_',alt1,'.mat'])

%% 1841-2010, subset 1961-2010

fyear = 1980:2010;
fyrall = 1841:2010;
fid = find(fyrall>=1980);

fmortD = fmortD(:,fid);
fmortF = fmortF(:,fid);
fmortP = fmortP(:,fid);

%% 1/2 degree
lats = unique([LatD; LatF; LatP]);
lons = unique([LonD; LonF; LonP]);

%%
nt = length(fid);
fmD = zeros(NID,nt);
fmF = zeros(NID,nt);
fmP = zeros(NID,nt);

for t=1:nt
    clear testD testF testP

    testD = griddata(LonD,LatD,fmortD(:,t),LON,LAT);
    fmD(:,t) = testD(WID);

    testF = griddata(LonF,LatF,fmortF(:,t),LON,LAT);
    fmF(:,t) = testF(WID);

    testP = griddata(LonP,LatP,fmortP(:,t),LON,LAT);
    fmP(:,t) = testP(WID);
end

%%
figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testF)
%caxis([0 1.2])
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testP)
%caxis([0 1.2])
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

figure
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,testD)
%caxis([0 1.2])
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% scale with Fmsy and temp
% fm = F/Fmsy, need to mult by Fmsy ~= 0.3
%tsc = (exp(0.063*(temp-10.0));

fmF = 0.3 * fmF .* (exp(0.063*(vmtp-10.0)));
fmP = 0.3 * fmP .* (exp(0.063*(vmtp-10.0)));
fmD = 0.3 * fmD .* (exp(0.063*(vmtb-10.0)));

%
fmD(isnan(fmD)) = 0.0;
fmF(isnan(fmF)) = 0.0;
fmP(isnan(fmP)) = 0.0;

fmD(fmD<0) = 0.0;
fmF(fmF<0) = 0.0;
fmP(fmP<0) = 0.0;

%% save
save([Cdir 'nemuro_ipsl_fmort_ID_annual_1980_2010_tempSc_',alt1,'.mat'],'year','WID',...
    'fmD','fmF','fmP');

Cdir2 ='/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
save([Cdir2 'nemuro_ipsl_fmort_ID_annual_1980_2010_tempSc_',alt1,'.mat'],'year','WID',...
    'fmD','fmF','fmP');



