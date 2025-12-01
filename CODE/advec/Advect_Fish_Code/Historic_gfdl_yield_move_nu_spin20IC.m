%%%%!! RUN SPINUP FOR ALL LOCATIONS
function Historic_gfdl_yield_move_nu_spin20IC()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;
param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_parameters_BCC(param);

%! Grids
%vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
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
param.adt = 12 * 60 * 60; %time step in seconds

%! How long to run the model
YEARS = 1988:2007;
nYEARS = length(YEARS);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = 'move_nu_v28_dt12h';
%opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
opath = '/project/Feisty/NC/Matlab_new_size/';
[fname,simname,sname] = sub_fname_hist_gfdl_move(param,opath,exper);

%! Storage variables
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

%! Initialize
load([sname '_' simname '.mat']); %Last month of spinup
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%! Dims of netcdf file
nt = 12 * nYEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% %%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_med_f = [fname,'_spin20IC_yield_med_f.nc'];
file_med_p = [fname,'_spin20IC_yield_med_p.nc'];
file_med_d = [fname,'_spin20IC_yield_med_d.nc'];
file_lrg_p = [fname,'_spin20IC_yield_lrg_p.nc'];
file_lrg_d = [fname,'_spin20IC_yield_lrg_d.nc'];

ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~15 minutes ... ']
xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidcatchMF    = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidcatchMP    = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidcatchMD    = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidcatchLP    = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidcatchLD    = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidLD,'time','double',time_dim);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model

addpath('matlab_functions');
MNT = 0;

%! Run model with no fishing
for YR = 1:nYEARS % years
    ti = num2str(YEARS(YR))
    load([vpath,'Data_ocean_cobalt_daily_',ti,'.mat'],'COBALT');
    load([vpath,'Vel100_esm2m_core_daily_',ti,'.mat'],'ESM');
    COBALT.U = ESM.U;
    COBALT.V = ESM.V;

    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_move_nu(DY,COBALT,GRD1,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param,neighborhood);

        %! Store
        S_Med_f(:,DY) = Med_f.caught;
        S_Med_p(:,DY) = Med_p.caught;
        S_Med_d(:,DY) = Med_d.caught;
        S_Lrg_p(:,DY) = Lrg_p.caught;
        S_Lrg_d(:,DY) = Lrg_d.caught;

    end %Days


    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker

        %! Put vars of netcdf file
        netcdf.putVar(ncidMF,vidcatchMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidcatchMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidcatchMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidcatchLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidcatchLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidTB,MNT-1,1,MNT);

    end %Monthly mean

end %Years

%! Close save
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);


end
