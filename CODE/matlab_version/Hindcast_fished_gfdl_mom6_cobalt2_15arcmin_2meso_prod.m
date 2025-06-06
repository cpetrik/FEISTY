%%%%!! RUN HISTORIC FOR ALL LOCATIONS
function Hindcast_fished_gfdl_mom6_cobalt2_15arcmin_2meso_prod()

% run on rockfish, so "/Volumes/petrik-lab/" now "/project/"

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
cpath='/project/Feisty/GCM_Data/MOM6-COBALTv2_reanalysis/';
load([cpath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat'],'GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:param.NX;

%! How long to run the model
YEARS = 1961:2010;
nYEARS = length(YEARS);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
opath = '/project/Feisty/NC/Matlab_new_size/';
[fname,sname,simname] = sub_fname_hindcast_gfdl_15arcmin_2meso(param,opath);

%! Storage variables
S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Sml_d = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

%! Initialize
load([sname '_' simname '.mat']); %Last month of spinup
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_prod_sml_f.nc'];
file_sml_p = [fname,'_prod_sml_p.nc'];
file_sml_d = [fname,'_prod_sml_d.nc'];
file_med_f = [fname,'_prod_med_f.nc'];
file_med_p = [fname,'_prod_med_p.nc'];
file_med_d = [fname,'_prod_med_d.nc'];
file_lrg_p = [fname,'_prod_lrg_p.nc'];
file_lrg_d = [fname,'_prod_lrg_d.nc'];

ncidSF = netcdf.create(file_sml_f,'NETCDF4');
ncidSP = netcdf.create(file_sml_p,'NETCDF4');
ncidSD = netcdf.create(file_sml_d,'NETCDF4');
ncidMF = netcdf.create(file_med_f,'NETCDF4');
ncidMP = netcdf.create(file_med_p,'NETCDF4');
ncidMD = netcdf.create(file_med_d,'NETCDF4');
ncidLP = netcdf.create(file_lrg_p,'NETCDF4');
ncidLD = netcdf.create(file_lrg_d,'NETCDF4');

%! Dims of netcdf file
nt = 12*nYEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~10 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
vidT        = netcdf.defVar(ncidLD,'time','double',time_dim);
netcdf.endDef(ncidLD);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:nYEARS % years
    %! Load a year's ESM data
    ti = num2str(YEARS(YR));
    ti
    load([cpath,'Data_gfdl_mom6_cobalt2_2meso_15arcmin_daily_',ti,'.mat'],'ESM');

    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
%         [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

        %! Store
        S_Sml_f(:,DY) = Sml_f.prod;
        S_Sml_p(:,DY) = Sml_p.prod;
        S_Sml_d(:,DY) = Sml_d.prod;
        S_Med_f(:,DY) = Med_f.prod;
        S_Med_p(:,DY) = Med_p.prod;
        S_Med_d(:,DY) = Med_d.prod;
        S_Lrg_p(:,DY) = Lrg_p.prod;
        S_Lrg_d(:,DY) = Lrg_d.prod;

    end %Days

    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker

        %! Put vars of netcdf file
        netcdf.putVar(ncidLD,vidT,MNT-1,1,MNT);

        netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(S_Sml_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(S_Sml_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidprodSD,[0 MNT-1],[NX 1],mean(S_Sml_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidprodMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidprodLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));


    end %Monthly mean

end %Years

%%
%! Close save
netcdf.close(ncidSF);
netcdf.close(ncidSP);
netcdf.close(ncidSD);
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
