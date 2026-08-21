%%%%!! RUN SPINUP FOR ALL LOCATIONS
function Historic_fished_yield_gfdl_NEP_proj()

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
%vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/MOM6-NEP10/';
vpath = '/project/Feisty/GCM_Data/MOM6-NEP10/';

%1-D
load([vpath 'Data_grid_mom6_nep10.mat'],'GRD');
GRD1 = GRD;
%clear GRD

% %2-D
% load([vpath 'Data_hindcast_grid_cp2D.mat'],'GRD')
% GRD2 = GRD;
% clear GRD
% 
% %grid params
% [ni,nj] = size(GRD2.mask);
% param.ni = ni;
% param.nj = nj;
% param.dx = GRD2.dxtn;
% param.dy = GRD2.dyte;
% param.mask = GRD2.mask;

param.NX = length(GRD1.ID);
param.ID = 1:param.NX;
NX = length(GRD1.ID);
ID = 1:param.NX;

%! Advection/Movement time step
% param.adt = 24 * 60 * 60; %time step in seconds

%! How long to run the model
YEARS = 1993:2024;
nYEARS = length(YEARS);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = '1993_no_move';
%opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
opath = '/project/Feisty/NC/Matlab_new_size/';
[fname,simname,sname] = sub_fname_hist_gfdl_nep(param,opath,exper);

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
file_med_f = [fname,'_yield_med_f.nc'];
file_med_p = [fname,'_yield_med_p.nc'];
file_med_d = [fname,'_yield_med_d.nc'];
file_lrg_p = [fname,'_yield_lrg_p.nc'];
file_lrg_d = [fname,'_yield_lrg_d.nc'];

ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~10 minutes ... ']
xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidyieldMF    = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidyieldMP    = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidyieldMD    = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidyieldLP    = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidyieldLD  = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
vidT       = netcdf.defVar(ncidLD,'time','double',time_dim);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model

MNT = 0;
%! Run model with no fishing
for YR = 1:nYEARS % years
    ti = num2str(YEARS(YR))
    load([vpath,'Data_mom6_nep10_daily_',ti,'.mat'],'ESM');
    
	% COBALT.U = uh;
    % COBALT.V = vh;

    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,ESM,GRD1,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

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
        netcdf.putVar(ncidMF,vidyieldMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidyieldMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidyieldMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidyieldLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidyieldLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidT,MNT-1,1,MNT);

    end %Monthly mean

end %Years

%! Close save
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);



end
