%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Hist_nemuro_obsfish_pp()

esms = {'IPSL','GFDL','HAD'};
esms2 = {'ipsl','gfdl','hadley'};
for mod = 1:3
    vers = esms{mod};
    esm = esms2{mod};

    if mod==1
        gcpath = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/IPSLdown/';
    elseif mod==2
        gcpath = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/GFDLdown/';
    elseif mod==3
        gcpath = '/Volumes/petrik-lab/Feisty/GCM_Data/NEMURO/HADdown/';
    end
    

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
load([gcpath 'nemuro_',esm,'_fmort_ID_annual_1980_2010_tempSc_assessment.mat'],...
    'fmD','fmF','fmP');

% Set fishing rate as 1st year for fname
param.frate = nan;
param.frateF = fmF(:,1);
param.frateP = fmP(:,1);
param.frateD = fmD(:,1);
param.dfrateF = param.frateF/365.0;
param.dfrateP = param.frateP/365.0;
param.dfrateD = param.frateD/365.0;

%! Make core parameters/constants (global)
param = make_parameters(param);

%! How long to run the model
modyrs = 1980:2010;
YEARS = length(modyrs);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Grid (choose where and when to run the model)
load([gcpath 'Data_grid_nemuro_',esm,'.mat'],'GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = length(GRD.Z);
ID = 1:NX;

%! Create a directory for output
[fname,simname,outdir] = fname_hist(param,vers);

%! Storage variables
S_Med_f_pp = zeros(NX,DAYS);
S_Lrg_p_pp = zeros(NX,DAYS);
S_Lrg_d_pp = zeros(NX,DAYS);

%! Initialize
%!From a previous run
load([outdir 'Last_mo_Spinup_',vers,'_All_fishobs_',simname,'.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_med_f = [fname,'_pp_med_f.nc'];
file_lrg_p = [fname,'_pp_lrg_p.nc'];
file_lrg_d = [fname,'_pp_lrg_d.nc'];

ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidppMF    = netcdf.defVar(ncidMF,'popprod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidppLP    = netcdf.defVar(ncidLP,'popprod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidppLD    = netcdf.defVar(ncidLD,'popprod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Run model with no fishing
MNT=0;
for YR = 1:YEARS % years
    
    MY = num2str(modyrs(YR))
    load([gcpath 'Data_nemuro_',esm,'_',MY,'.mat'],'ESM');

    param.frateF = fmF(:,YR);
    param.frateP = fmP(:,YR);
    param.frateD = fmD(:,YR);
    param.dfrateF = param.frateF/365.0;
    param.dfrateP = param.frateP/365.0;
    param.dfrateD = param.frateD/365.0;
    
    for DY = 1:DAYS % days


        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

        S_Med_f_pp(:,DY) = Med_f.pp;
        S_Lrg_p_pp(:,DY) = Lrg_p.pp;
        S_Lrg_d_pp(:,DY) = Lrg_d.pp;

    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker

        %! Put vars of netcdf file
        netcdf.putVar(ncidMF,vidppMF,[0 MNT-1],[NX 1],mean(S_Med_f_pp(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidppLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_pp(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidppLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_pp(:,a(i):b(i)),2));

    end %Monthly mean

end %Years

%! Close save
netcdf.close(ncidMF);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
