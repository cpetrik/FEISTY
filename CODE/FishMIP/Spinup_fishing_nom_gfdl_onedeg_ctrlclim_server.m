%%%%!! RUN SPINUP GLOBALLY
% GFDL reanalysis, ctrlclim, 1 degree
function Spinup_fishing_nom_gfdl_onedeg_ctrlclim_server()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
vers = 'nominal';
% Nominal effort method of estimating
load(['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',vers,'/',...
    'gfdl-mom6-cobalt2_ctrlclim_onedeg_fmort_ID_annual_1841_1960_tempSc_',vers,'.mat'],...
    'fmD','fmP','fmF');

% Set fishing rate as 1st year 1841
param.frate = nan;
param.frateF = fmF(:,1);
param.frateP = fmP(:,1);
param.frateD = fmD(:,1);
param.dfrateF = param.frateF/365.0;
param.dfrateP = param.frateP/365.0;
param.dfrateD = param.frateD/365.0;

%! Make core parameters/constants
param = make_parameters_1meso(param);

%! Grid
load('/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat','GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:param.NX;

%! How long to run the model
CYCLES = 10;
YEARS = 1961:1980;
nYEARS = length(YEARS);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
[fname,simname] = sub_fname_spin_gfdl_onedeg_vers_server(param,vers);

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

C_Med_f = zeros(NX,DAYS);
C_Med_p = zeros(NX,DAYS);
C_Med_d = zeros(NX,DAYS);
C_Lrg_p = zeros(NX,DAYS);
C_Lrg_d = zeros(NX,DAYS);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(param.ID,DAYS);
%ENVR = sub_init_env_empHP(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_empHP_sml_f.nc'];
file_sml_p = [fname,'_empHP_sml_p.nc'];
file_sml_d = [fname,'_empHP_sml_d.nc'];
file_med_f = [fname,'_empHP_med_f.nc'];
file_med_p = [fname,'_empHP_med_p.nc'];
file_med_d = [fname,'_empHP_med_d.nc'];
file_lrg_p = [fname,'_empHP_lrg_p.nc'];
file_lrg_d = [fname,'_empHP_lrg_d.nc'];
file_bent  = [fname,'_empHP_bent.nc'];
% file_mzoo  = [fname,'_mzoo.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');
ncidB  = netcdf.create(file_bent,'NC_WRITE');
% ncidMZ = netcdf.create(file_mzoo,'NC_WRITE');

%! Dims of netcdf file
nt = 12 * nYEARS * CYCLES;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~10 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
vidcatchMF  = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
vidcatchMP  = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
vidcatchMD  = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
vidcatchLP  = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
vidcatchLD  = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

xy_dim     = netcdf.defDim(ncidB,'nid',NX);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for c = 1:CYCLES
    for YR = 1:nYEARS % years
        %! Load a year's CESM data
        ti = num2str(YEARS(YR));
        [ti,' , ', num2str(c)]
        load(['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/',...
            'Data_gfdl_mom6_cobalt2_ctrlclim_onedeg_daily_',ti,'.mat'],'ESM');

        for DAY = 1:DAYS % days

            %%%! Future time step
            DY = int64(ceil(DAY));
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                sub_futbio_1meso_empHPloss_obsfished(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
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

            C_Med_f(:,DY) = Med_f.caught;
            C_Med_p(:,DY) = Med_p.caught;
            C_Med_d(:,DY) = Med_d.caught;
            C_Lrg_p(:,DY) = Lrg_p.caught;
            C_Lrg_d(:,DY) = Lrg_d.caught;

        end %Days

        %! Calculate monthly means and save
        %if (c==CYCLES) % save last 20 yrs
        aa = (cumsum(MNTH)+1);
        a = [1,aa(1:end-1)]; % start of the month
        b = cumsum(MNTH); % end of the month
        for i = 1:12
            MNT = MNT+1; % Update monthly ticker

            %! Put vars of netcdf file
            netcdf.putVar(ncidB,vidbioB,[0 MNT-1],[NX 1],mean(S_Bent_bio(:,a(i):b(i)),2));
            netcdf.putVar(ncidB,vidTB,MNT-1,1,MNT);
          
            netcdf.putVar(ncidSF,vidbioSF,[0 MNT-1],[NX 1],mean(S_Sml_f(:,a(i):b(i)),2));
            netcdf.putVar(ncidSP,vidbioSP,[0 MNT-1],[NX 1],mean(S_Sml_p(:,a(i):b(i)),2));
            netcdf.putVar(ncidSD,vidbioSD,[0 MNT-1],[NX 1],mean(S_Sml_d(:,a(i):b(i)),2));
            netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
            netcdf.putVar(ncidMP,vidbioMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
            netcdf.putVar(ncidMD,vidbioMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
            netcdf.putVar(ncidLP,vidbioLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
            netcdf.putVar(ncidLD,vidbioLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));

            netcdf.putVar(ncidMF,vidcatchMF,[0 MNT-1],[NX 1],mean(C_Med_f(:,a(i):b(i)),2));
            netcdf.putVar(ncidMP,vidcatchMP,[0 MNT-1],[NX 1],mean(C_Med_p(:,a(i):b(i)),2));
            netcdf.putVar(ncidMD,vidcatchMD,[0 MNT-1],[NX 1],mean(C_Med_d(:,a(i):b(i)),2));
            netcdf.putVar(ncidLP,vidcatchLP,[0 MNT-1],[NX 1],mean(C_Lrg_p(:,a(i):b(i)),2));
            netcdf.putVar(ncidLD,vidcatchLD,[0 MNT-1],[NX 1],mean(C_Lrg_d(:,a(i):b(i)),2));

        end %Monthly mean
        %end

    end %Years
end %Cycles

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
netcdf.close(ncidB);

end