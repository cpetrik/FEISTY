%%%%!! RUN PRE-INDUSTRIAL FOR ALL LOCATIONS
function PI_fishing_empHP_gfdl_mom6_cobalt2_15arcmin_ctrlclim_server_MF()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
% load(['/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/',...
%     'gfdl-mom6-cobalt2_ctrlclim_15arcmin_fmort_ID_annual_1841_1960_tempSc.mat'],...
%     'fmD','fmF','fmP');
load(['/Users/cpetrik/Dropbox/Princeton/FEISTY_other/fishing_ms_ideas/fishing_effort_impl/grid_mortality_guilds_v1/',...
    'gfdl-mom6-cobalt2_ctrlclim_15arcmin_fmort_ID_annual_1841_1960_tempSc.mat'],...
    'fmD','fmF','fmP');
% Set fishing rate as 1st year for fname
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
%load('/Volumes/MIP/Fish-MIP/Phase3/QuarterDeg/Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat','GRD');
load('/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_15arcmin.mat','GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:param.NX;

%! How long to run the model
CYCLES = 6;
YEARS = 1961:1980;
nYEARS = length(YEARS);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
[fname,simname,outdir] = sub_fname_pi_gfdl_15arcmin_ctrl_server(param);

%! Storage variables
S_Med_f = zeros(NX,DAYS);

%! Initialize
load([fname '_Last_mo_' simname '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

for c = 1:CYCLES
    %%%%%%%%%%%%%%% Setup NetCDF save
    %! Setup netcdf path to store to
    file_med_f  = [fname,'_empHP_med_f_cycle',num2str(c),'.nc'];

    ncidMF  = netcdf.create(file_med_f,'NC_WRITE');

    %! Dims of netcdf file
    nt = 12 * nYEARS;
    netcdf.setDefaultFormat('NC_FORMAT_64BIT');

    %% ! Def vars of netcdf file
    ['Defining netcdfs, takes ~5 minutes ... ']
    xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
    time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
    vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
    netcdf.endDef(ncidMF);


    %% %%%%%%%%%%%%%%%%%%%% Run the Model
    MNT = 0;
    FYR = 0;
    %! Run model with fishing

    for YR = 1:nYEARS % years
        %! Load a year's ESM data
        ti = num2str(YEARS(YR));
        [ti,' , ', num2str(c)]
        load(['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/QuarterDeg/',...
            'Data_gfdl_mom6_cobalt2_ctrlclim_15arcmin_daily_',ti,'.mat'],'ESM');

        FYR = FYR+1; %Climate cycles, but fishing varies by yr
        param.frateF = fmF(:,FYR);
        param.frateP = fmP(:,FYR);
        param.frateD = fmD(:,FYR);
        param.dfrateF = param.frateF/365.0;
        param.dfrateP = param.frateP/365.0;
        param.dfrateD = param.frateD/365.0;

        for DAY = 1:param.DT:DAYS % days

            %%%! Future time step
            DY = int64(ceil(DAY));
            %         [num2str(YR),' , ', num2str(mod(DY,365))]
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                sub_futbio_1meso_empHPloss_obsfished(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
                Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

            %! Store
            S_Med_f(:,DY) = Med_f.bio;

        end %Days

        %! Calculate monthly means and save
        aa = (cumsum(MNTH)+1);
        a = [1,aa(1:end-1)]; % start of the month
        b = cumsum(MNTH); % end of the month
        for i = 1:12
            MNT = MNT+1; % Update monthly ticker

            %! Put vars of netcdf file
            netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));

        end %Monthly mean

    end %Years


    %%
    %! Close save
    netcdf.close(ncidMF);

end %cycles

end
