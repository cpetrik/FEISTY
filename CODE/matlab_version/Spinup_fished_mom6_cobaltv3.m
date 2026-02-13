%%%%!! RUN SPINUP FOR ALL LOCATIONS
function Spinup_fished_mom6_cobaltv3()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;
param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_parameters(param);

%! Grid
% cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
% gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
cpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = cpath;
load([cpath 'Data_grid_OM4_05_COBALTv3.mat'],'GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:param.NX;

%! How long to run the model
YEARS = 200;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = 'Spinup1990';
%opath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/';
opath = '/project/Feisty/NC/Matlab_new_size/';
[fname,simname,outdir] = sub_fname_spin_mom6_cobaltv2(param,opath,exper);

%! Storage variables
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

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);

%! Dims of netcdf file
nt = (YEARS);
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% %%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_sml_f.nc'];
file_sml_p = [fname,'_sml_p.nc'];
file_sml_d = [fname,'_sml_d.nc'];
file_med_f = [fname,'_med_f.nc'];
file_med_p = [fname,'_med_p.nc'];
file_med_d = [fname,'_med_d.nc'];
file_lrg_p = [fname,'_lrg_p.nc'];
file_lrg_d = [fname,'_lrg_d.nc'];
file_bent  = [fname,'_bent.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');
ncidB  = netcdf.create(file_bent,'NC_WRITE');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~15 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
% vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
% vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,time_dim]);
% vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
% vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
% vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
% vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
% vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
% vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

xy_dim     = netcdf.defDim(ncidB,'nid',NX);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);

%% %%%%%%%%%%%%%%%%%%%% Run the Model

load([gpath,'ocean_cobalt_feisty_forcing_2dint.1990.mat'],'ESM');

MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    ti = num2str(YR)

    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

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


    %! Calculate annual means and save
    %! Put vars of netcdf file
    netcdf.putVar(ncidB,vidbioB,[0 YR-1],[NX 1],mean(S_Bent_bio,2));
    netcdf.putVar(ncidB,vidTB,YR-1,1,YR);

    netcdf.putVar(ncidSF,vidbioSF,[0 YR-1],[NX 1],mean(S_Sml_f,2));
    netcdf.putVar(ncidSP,vidbioSP,[0 YR-1],[NX 1],mean(S_Sml_p,2));
    netcdf.putVar(ncidSD,vidbioSD,[0 YR-1],[NX 1],mean(S_Sml_d,2));
    netcdf.putVar(ncidMF,vidbioMF,[0 YR-1],[NX 1],mean(S_Med_f,2));
    netcdf.putVar(ncidMP,vidbioMP,[0 YR-1],[NX 1],mean(S_Med_p,2));
    netcdf.putVar(ncidMD,vidbioMD,[0 YR-1],[NX 1],mean(S_Med_d,2));
    netcdf.putVar(ncidLP,vidbioLP,[0 YR-1],[NX 1],mean(S_Lrg_p,2));
    netcdf.putVar(ncidLD,vidbioLD,[0 YR-1],[NX 1],mean(S_Lrg_d,2));


end %Years

%! Close save
netcdf.close(ncidSF);
netcdf.close(ncidSP);
netcdf.close(ncidSD);
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);
netcdf.close(ncidB);



end
