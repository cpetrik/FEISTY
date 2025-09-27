%%%%!! RUN SPINUP FOR ALL LOCATIONS
function Spinup_fished_gfdl_move_nu_Y1()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;
param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_parameters(param);

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
YEARS = 1; %30;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
%opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
opath = '/project/Feisty/NC/Matlab_new_size/';
exper = 'Spinup1988_move_nu_v24_dt12h';
[fname,simname,sname] = sub_fname_spin_move_core(param,opath,exper);

%! Storage variables
S_Bent_bio = zeros(NX,DAYS);

S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Sml_d = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

%! Initialize
%[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);

% Last month of spinup without movement
load([sname '_' simname '.mat']); 
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%! Dims of netcdf file
nt = 12 * YEARS;


%% %%%%%%%%%%%%%%%%%%%% Run the Model

addpath('matlab_functions');

load([vpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');
load([vpath,'Vel100_esm2m_core_daily_1988.mat'],'ESM');
COBALT.U = ESM.U;
COBALT.V = ESM.V;

MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    ti = num2str(YR)

    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_move_nu(DY,COBALT,GRD1,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param,neighborhood);

        %! Store
        S_Bent_bio(:,DY) = BENT.mass;

        S_Sml_f(:,DY) = Sml_f.bio;
        S_Sml_p(:,DY) = Sml_p.bio;
        S_Sml_d(:,DY) = Sml_d.bio;
        S_Med_f(:,DY) = Med_f.bio;
        S_Med_p(:,DY) = Med_p.bio;
        S_Med_d(:,DY) = Med_d.bio;
        S_Lrg_p(:,DY) = Lrg_p.bio;
        S_Lrg_d(:,DY) = Lrg_d.bio;

    end %Days

end %Years

save([fname,'_Y1.mat'],'S_Bent_bio','S_Sml_f','S_Sml_p','S_Sml_d','S_Med_f',...
    'S_Med_p','S_Med_d','S_Lrg_p','S_Lrg_d','GRD1','GRD2','exper');


end
