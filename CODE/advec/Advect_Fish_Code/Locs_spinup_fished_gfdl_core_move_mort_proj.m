%%%%!! RUN SPINUP FOR ONE LOCATION
function Locs_spinup_fished_gfdl_core_move_mort_proj()

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

%Single locations
load([vpath 'core_grid_360x200_id_locs_area_dep.mat'],'ids','abbrev');
names = abbrev;

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
param.adt = 6 * 60 * 60; %time step in seconds

%! How long to run the model
YEARS = 50; 
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = 'Locs_spinup1988_move_mort_v21_dt6h';
%opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
opath = '/project/Feisty/NC/Matlab_new_size/';
[fname,simname,sname] = sub_fname_spin_move_core(param,opath,exper);

%! Initialize
%[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);

% Last month of spinup without movement
load([sname '_' simname '.mat']); 
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%! Storage
nid = length(ids);
Spinup_Sml_f = NaN*ones(DAYS,25,nid);
Spinup_Sml_p = NaN*ones(DAYS,25,nid);
Spinup_Sml_d = NaN*ones(DAYS,25,nid);
Spinup_Med_f = NaN*ones(DAYS,25,nid);
Spinup_Med_p = NaN*ones(DAYS,25,nid);
Spinup_Med_d = NaN*ones(DAYS,25,nid);
Spinup_Lrg_p = NaN*ones(DAYS,25,nid);
Spinup_Lrg_d = NaN*ones(DAYS,25,nid);
Spinup_Cobalt = NaN*ones(DAYS,5,nid);

S_Sml_f = NaN*ones(12*YEARS,25,nid);
S_Sml_p = NaN*ones(12*YEARS,25,nid);
S_Sml_d = NaN*ones(12*YEARS,25,nid);
S_Med_f = NaN*ones(12*YEARS,25,nid);
S_Med_p = NaN*ones(12*YEARS,25,nid);
S_Med_d = NaN*ones(12*YEARS,25,nid);
S_Lrg_p = NaN*ones(12*YEARS,25,nid);
S_Lrg_d = NaN*ones(12*YEARS,25,nid);
S_Cobalt = NaN*ones(12*YEARS,5,nid);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
addpath('matlab_functions');

load([vpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');
load([vpath,'Vel100_esm2m_core_daily_1988.mat'],'ESM');
COBALT.U = ESM.U;
COBALT.V = ESM.V;

MNT=0;
for YR = 1:YEARS % years
    num2str(YR)

    for DAY = 1:DAYS % days
        
        %%%! ticker
        DY = int64(ceil(DAY));
        
        %%%! Future time step
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_move_mort(DY,COBALT,GRD1,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param,neighborhood);
        
        %! Store last year
        %if (YR==YEARS)
        Spinup_Sml_f(DY,:,:) = [Sml_f.bio(ids) Sml_f.enc_f(ids) Sml_f.enc_p(ids) Sml_f.enc_d(ids) Sml_f.enc_zm(ids) ...
            Sml_f.enc_zl(ids) Sml_f.enc_be(ids) Sml_f.con_f(ids) Sml_f.con_p(ids) Sml_f.con_d(ids) Sml_f.con_zm(ids) ...
            Sml_f.con_zl(ids) Sml_f.con_be(ids) Sml_f.I(ids) Sml_f.nu(ids) Sml_f.gamma(ids) Sml_f.die(ids) Sml_f.rep(ids) ...
            Sml_f.rec(ids) Sml_f.clev(ids) Sml_f.prod(ids) Sml_f.pred(ids) Sml_f.nmort(ids) Sml_f.met(ids) Sml_f.caught(ids)]';
        Spinup_Sml_p(DY,:,:) = [Sml_p.bio(ids) Sml_p.enc_f(ids) Sml_p.enc_p(ids) Sml_p.enc_d(ids) Sml_p.enc_zm(ids) Sml_p.enc_zl(ids) Sml_p.enc_be(ids) Sml_p.con_f(ids) Sml_p.con_p(ids) Sml_p.con_d(ids) Sml_p.con_zm(ids) Sml_p.con_zl(ids) Sml_p.con_be(ids) Sml_p.I(ids) Sml_p.nu(ids) Sml_p.gamma(ids) Sml_p.die(ids) Sml_p.rep(ids) Sml_p.rec(ids) Sml_p.clev(ids) Sml_p.prod(ids) Sml_p.pred(ids) Sml_p.nmort(ids) Sml_p.met(ids) Sml_p.caught(ids)]';
        Spinup_Sml_d(DY,:,:) = [Sml_d.bio(ids) Sml_d.enc_f(ids) Sml_d.enc_p(ids) Sml_d.enc_d(ids) Sml_d.enc_zm(ids) Sml_d.enc_zl(ids) Sml_d.enc_be(ids) Sml_d.con_f(ids) Sml_d.con_p(ids) Sml_d.con_d(ids) Sml_d.con_zm(ids) Sml_d.con_zl(ids) Sml_d.con_be(ids) Sml_d.I(ids) Sml_d.nu(ids) Sml_d.gamma(ids) Sml_d.die(ids) Sml_d.rep(ids) Sml_d.rec(ids) Sml_d.clev(ids) Sml_d.prod(ids) Sml_d.pred(ids) Sml_d.nmort(ids) Sml_d.met(ids) Sml_d.caught(ids)]';
        Spinup_Med_f(DY,:,:) = [Med_f.bio(ids) Med_f.enc_f(ids) Med_f.enc_p(ids) Med_f.enc_d(ids) Med_f.enc_zm(ids) Med_f.enc_zl(ids) Med_f.enc_be(ids) Med_f.con_f(ids) Med_f.con_p(ids) Med_f.con_d(ids) Med_f.con_zm(ids) Med_f.con_zl(ids) Med_f.con_be(ids) Med_f.I(ids) Med_f.nu(ids) Med_f.gamma(ids) Med_f.die(ids) Med_f.rep(ids) Med_f.rec(ids) Med_f.clev(ids) Med_f.prod(ids) Med_f.pred(ids) Med_f.nmort(ids) Med_f.met(ids) Med_f.caught(ids)]';
        Spinup_Med_p(DY,:,:) = [Med_p.bio(ids) Med_p.enc_f(ids) Med_p.enc_p(ids) Med_p.enc_d(ids) Med_p.enc_zm(ids) Med_p.enc_zl(ids) Med_p.enc_be(ids) Med_p.con_f(ids) Med_p.con_p(ids) Med_p.con_d(ids) Med_p.con_zm(ids) Med_p.con_zl(ids) Med_p.con_be(ids) Med_p.I(ids) Med_p.nu(ids) Med_p.gamma(ids) Med_p.die(ids) Med_p.rep(ids) Med_p.rec(ids) Med_p.clev(ids) Med_p.prod(ids) Med_p.pred(ids) Med_p.nmort(ids) Med_p.met(ids) Med_p.caught(ids)]';
        Spinup_Med_d(DY,:,:) = [Med_d.bio(ids) Med_d.enc_f(ids) Med_d.enc_p(ids) Med_d.enc_d(ids) Med_d.enc_zm(ids) Med_d.enc_zl(ids) Med_d.enc_be(ids) Med_d.con_f(ids) Med_d.con_p(ids) Med_d.con_d(ids) Med_d.con_zm(ids) Med_d.con_zl(ids) Med_d.con_be(ids) Med_d.I(ids) Med_d.nu(ids) Med_d.gamma(ids) Med_d.die(ids) Med_d.rep(ids) Med_d.rec(ids) Med_d.clev(ids) Med_d.prod(ids) Med_d.pred(ids) Med_d.nmort(ids) Med_d.met(ids) Med_d.caught(ids)]';
        Spinup_Lrg_p(DY,:,:) = [Lrg_p.bio(ids) Lrg_p.enc_f(ids) Lrg_p.enc_p(ids) Lrg_p.enc_d(ids) Lrg_p.enc_zm(ids) Lrg_p.enc_zl(ids) Lrg_p.enc_be(ids) Lrg_p.con_f(ids) Lrg_p.con_p(ids) Lrg_p.con_d(ids) Lrg_p.con_zm(ids) Lrg_p.con_zl(ids) Lrg_p.con_be(ids) Lrg_p.I(ids) Lrg_p.nu(ids) Lrg_p.gamma(ids) Lrg_p.die(ids) Lrg_p.rep(ids) Lrg_p.rec(ids) Lrg_p.clev(ids) Lrg_p.prod(ids) Lrg_p.pred(ids) Lrg_p.nmort(ids) Lrg_p.met(ids) Lrg_p.caught(ids)]';
        Spinup_Lrg_d(DY,:,:) = [Lrg_d.bio(ids) Lrg_d.enc_f(ids) Lrg_d.enc_p(ids) Lrg_d.enc_d(ids) Lrg_d.enc_zm(ids) Lrg_d.enc_zl(ids) Lrg_d.enc_be(ids) Lrg_d.con_f(ids) Lrg_d.con_p(ids) Lrg_d.con_d(ids) Lrg_d.con_zm(ids) Lrg_d.con_zl(ids) Lrg_d.con_be(ids) Lrg_d.I(ids) Lrg_d.nu(ids) Lrg_d.gamma(ids) Lrg_d.die(ids) Lrg_d.rep(ids) Lrg_d.rec(ids) Lrg_d.clev(ids) Lrg_d.prod(ids) Lrg_d.pred(ids) Lrg_d.nmort(ids) Lrg_d.met(ids) Lrg_d.caught(ids)]';
        Spinup_Cobalt(DY,:,:) = [BENT.mass(ids) BENT.pred(ids) ENVR.fZm(ids) ENVR.fZl(ids) ENVR.fB(ids)]';
        %end
        
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        S_Cobalt(MNT,:,:) = mean(Spinup_Cobalt(a(i):b(i),:,:),1);
        S_Sml_f(MNT,:,:) = mean(Spinup_Sml_f(a(i):b(i),:,:),1);
        S_Sml_p(MNT,:,:) = mean(Spinup_Sml_p(a(i):b(i),:,:),1);
        S_Sml_d(MNT,:,:) = mean(Spinup_Sml_d(a(i):b(i),:,:),1);
        S_Med_f(MNT,:,:) = mean(Spinup_Med_f(a(i):b(i),:,:),1);
        S_Med_p(MNT,:,:) = mean(Spinup_Med_p(a(i):b(i),:,:),1);
        S_Med_d(MNT,:,:) = mean(Spinup_Med_d(a(i):b(i),:,:),1);
        S_Lrg_p(MNT,:,:) = mean(Spinup_Lrg_p(a(i):b(i),:,:),1);
        S_Lrg_d(MNT,:,:) = mean(Spinup_Lrg_d(a(i):b(i),:,:),1);
    end
    
end %Years

%%% Save
save([fname '.mat'],...
'S_Sml_f','S_Sml_p','S_Sml_d','S_Med_f','S_Med_p','S_Med_d',...
'S_Lrg_p','S_Lrg_d','S_Cobalt')

end
