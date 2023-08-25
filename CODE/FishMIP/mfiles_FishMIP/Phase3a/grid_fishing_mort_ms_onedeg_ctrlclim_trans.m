% Put fishing mortality onto grid
% 3 versions for fishing effort -> mort ms

clear
close all

%% grid info
% ctrlclim one degree
Cdir = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

% Depth, lat, lon, area, grid cell with seafloor
load([Cdir 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);
load([Cdir 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat']);

[ni,nj] = size(LON);

WID = GRD.ID;
NID = GRD.N;

clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;

%% temp scaling
load([Cdir 'gfdl-mom6-cobalt2_ctrlclim_temp100_onedeg_global_monthly_1961_2010.mat'],'temp_100');
load([Cdir 'gfdl-mom6-cobalt2_ctrlclim_tob_onedeg_global_monthly_1961_2010.mat'],'tob','yr');

% take mean temp over hist time period = 1961-2010
mtp = nanmean(temp_100,3);
mtb = nanmean(tob,3);

vmtp = mtp(WID);
vmtb = mtb(WID);

%% time
year = 1841:1960;
yrall = 1841:2010;
yid = find(yrall<=1960);

%% load f/fmsy
alt = {'assessment','effective','nominal'};

% loop through
for i=1:length(alt)
    alt1 = alt{i};
    spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
    fpath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
        'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
    load([fpath 'grid_mortality_all_',alt1,'.mat'])

    %% 1841-2010, subset 1841-1960
    fmortD = fmortD(:,yid);
    fmortF = fmortF(:,yid);
    fmortP = fmortP(:,yid);

    %% 1/2 degree
    lats = unique([LatD; LatF; LatP]);
    lons = unique([LonD; LonF; LonP]);

    %%
    nt = length(yid);
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
    caxis([0 0.6])
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    figure
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,testP)
    caxis([0 0.6])
    h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

    figure
    axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
        'Grid','off','FLineWidth',1)
    surfm(LAT,LON,testD)
    caxis([0 0.6])
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
    tpath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
    save([tpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc_',alt1,'.mat'],'year','WID',...
        'fmD','fmF','fmP');
    save([fpath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc_',alt1,'.mat'],'year','WID',...
        'fmD','fmF','fmP');
    save([spath 'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc_',alt1,'.mat'],'year','WID',...
        'fmD','fmF','fmP');

end


