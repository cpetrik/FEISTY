

vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
param.frate = 0.3;
param.dfrate = param.frate/365.0;
param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_parameters(param);

%! Grids
vpath = '/project/Feisty/GCM_Data/CORE-forced/';

%1-D
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
GRD1 = GRD;
clear GRD

%2-D
load([vpath 'Data_hindcast_grid_cp2D.mat'],'GRD')
GRD2 = GRD;
clear GRD

%Grid cell neighbors
load([vpath 'all_neighbors_core_grid_360x200.mat'],'neighborhood')

%grid params
[ni,nj] = size(GRD2.mask);
param.ni = ni;
param.nj = nj;
param.dx = GRD2.dxtn;
param.dy = GRD2.dyte;
param.mask = GRD2.mask;
param.area = GRD2.area;

param.NX = length(GRD1.Z);
param.ID = 1:param.NX;
NX = length(GRD1.Z);
ID = 1:param.NX;

%! Advection/Movement time step
param.adt = 24 * 60 * 60; %time step in seconds

%! How long to run the model
YEARS = 1; %50;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! How long to run the model
YEARS = 1; %50;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
exper = 'Spinup1988_move_prey_v3';
%[fname,simname,sname] = sub_fname_spin_gfdl_core(param,opath,exper);
[fname,simname,sname] = sub_fname_spin_move_core(param,opath,exper);
S_Bent_bio = zeros(NX,DAYS);
% S_Mzoo_frac = zeros(NX,DAYS);

S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Sml_d = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);
load([sname '_' simname '.mat']); 
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

addpath('matlab_functions');

load([vpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');
load([vpath,'Vel200_feb152013_run25_ocean_1988.mat'],'uh','vh');
COBALT.U = uh;
COBALT.V = vh;

Sf = Sml_f;
Sp = Sml_p;
Sd = Sml_d;
Mf = Med_f;
Mp = Med_p;
Md = Med_d;
Lp = Lrg_p;
Ld = Lrg_d;

YR=1;
DY=1;
ESM=COBALT;
GRD=GRD1;
neighbor=neighborhood;


preySf = sub_1Dto2D(GRD,Sf.I,param);
preySp = sub_1Dto2D(GRD,Sp.I,param);
preySd = sub_1Dto2D(GRD,Sd.I,param);
preyMf = sub_1Dto2D(GRD,Mf.I,param);
preyMp = sub_1Dto2D(GRD,Mp.I,param);
preyMd = sub_1Dto2D(GRD,Md.I,param);
preyLp = sub_1Dto2D(GRD,Lp.I,param);
preyLd = sub_1Dto2D(GRD,Ld.I,param);

bioSf = sub_1Dto2D(GRD,Sf.bio,param);
bioSp = sub_1Dto2D(GRD,Sp.bio,param);
bioSd = sub_1Dto2D(GRD,Sd.bio,param);
bioMf = sub_1Dto2D(GRD,Mf.bio,param);
bioMp = sub_1Dto2D(GRD,Mp.bio,param);
bioMd = sub_1Dto2D(GRD,Md.bio,param);
bioLp = sub_1Dto2D(GRD,Lp.bio,param);
bioLd = sub_1Dto2D(GRD,Ld.bio,param);


[m,n] = ind2sub([ni,nj],GRD.ID(13))

preyMf(m,n)

bioMf(m,n)


u200 = ENVR.U;
v200 = ENVR.V;

current = nan*ones(param.ni,param.nj,2);
current(:,:,1) = u200; 
current(:,:,2) = v200;
 
bioMf2 = AdvectPredator(bioMf,preyMf,current,param.adt,param.dx,param.dy,neighbor,param.U_m,param.mask,param.area,param.nj,param.ni);
bioMf2(m,n)


load('/Volumes/petrik-lab/Feisty/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');

[ni,nj]=size(geolon_t);
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; 
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,bioMf)


figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,bioMf2)

figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,u200)

figure(4)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,v200)


