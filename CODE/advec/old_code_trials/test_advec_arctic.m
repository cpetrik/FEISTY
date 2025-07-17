% Test advection in MOM6 grid at Arctic crossover

clear all
close all

% Transports path
vpath = '/Volumes/GFDL/GCM_DATA/MOM6/Preindust/';

% Grid
load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid1D_MOM6_preindust.mat','GRD');
ID = GRD.ID;
clear GRD
load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid2D_MOM6_preindust.mat')

tmask = GRD.lmask;

%% number of water cells
NX = length(ID);

% grid size
[ni,nj] = size(tmask);
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

%% define a patch to advect
bio = zeros(ni,nj);
%File name to save

%seed Arctic
bio(125:135,(jed-1):jed) = 1.0e2;    
cname='Arc_test_dt1hr_mom6_preind_vel_surf';


bio = bio .* tmask;

% define time
YEARS = 1;
DAYS = 365;
tstep = 1; %time step in hours

biov = zeros(NX,DAYS*YEARS);
bio1 = NaN*ones(ni,nj,DAYS*YEARS);

%% do advec only
% 1ST ORDER
K=0;
n=0;
for Y=1:YEARS
    % Velocities
    load([vpath 'MOM6_00010101_velocities_daily_interp.mat'],'D_u1','D_v1');
    U = zeros(ni,nj);
    V = zeros(ni,nj);
    for DAY = 1:3%:DAYS
        [num2str(DAY) ',1']
        n=n+1;
        u = D_u1(:,DAY);
        v = D_v1(:,DAY);
        U(ID) = u; 
        V(ID) = v;
        bio = sub_mom6_advec1_diff_vel(GRD,bio,K,U,V,ni,nj,tstep);
        biov(:,n) = bio(ID);
        bio1(:,:,n) = bio;
    end
end

% Save
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_1st_2nd_4th_' cname '_UV1.mat'],...
%     'biov','bio');

%%
[nid,nd] = size(biov);

% % Grid with area instead of vectors of area
% bio1 = NaN*ones(ni,nj,nd);
% for d = 1:nd
%     biot = NaN*ones(ni,nj);
%     biot(ID) = (biov(:,d));
%     bio1(:,:,d) = biot;
% end

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% Land
surf_tmask = tmask(:,:,1);
lmask = surf_tmask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

% colors
cmB=cbrewer('seq','Blues',50);

%% Global flat
for t=1:3
    B1 = bio1(:,:,t);
    
    figure
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(GRD.LAT,GRD.LON,B1)
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar;
    caxis([50 100]);
    colormap(cmB)
    title(['Day ' num2str(t) ' Year 1'])
    %print('-dpng',[fpath 'advec_diff_test_1st_K100_' cname '_' num2str(t(n)) '.png'])
end


%% Arctic projection
for t=1:3
    B1 = bio1(:,:,t);
    
    figure
    m_proj('stereographic','lat',90,'long',30,'radius',20);
    m_pcolor(GRD.LON,GRD.LAT,B1);
    shading interp
    colorbar
    colormap(cmB)
    caxis([0 100])
    m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t) ' Year 1'])
    %print('-dpng',[fpath 'advec_diff_test_1st_K100_' cname '_arcticproj_' num2str(t(n)) '.png'])
    
end
