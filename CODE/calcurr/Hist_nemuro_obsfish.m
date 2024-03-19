%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Hist_nemuro_obsfish()

esms = {'IPSL','GFDL','HAD'};
esms2 = {'ipsl','gfdl','hadley'};
for mod = 3
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
S_Bent_bio = zeros(NX,DAYS);

S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Sml_d = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

S_Med_f_fish = zeros(NX,DAYS);
S_Med_p_fish = zeros(NX,DAYS);
S_Med_d_fish = zeros(NX,DAYS);
S_Lrg_p_fish = zeros(NX,DAYS);
S_Lrg_d_fish = zeros(NX,DAYS);

% S_Sml_f_rec = zeros(NX,DAYS);
% S_Sml_p_rec = zeros(NX,DAYS);
% S_Sml_d_rec = zeros(NX,DAYS);
% S_Med_f_rec = zeros(NX,DAYS);
% S_Med_p_rec = zeros(NX,DAYS);
% S_Med_d_rec = zeros(NX,DAYS);
% S_Lrg_p_rec = zeros(NX,DAYS);
% S_Lrg_d_rec = zeros(NX,DAYS);
% 
% S_Sml_f_con = zeros(NX,DAYS);
% S_Sml_p_con = zeros(NX,DAYS);
% S_Sml_d_con = zeros(NX,DAYS);
% S_Med_f_con = zeros(NX,DAYS);
% S_Med_p_con = zeros(NX,DAYS);
% S_Med_d_con = zeros(NX,DAYS);
% S_Lrg_p_con = zeros(NX,DAYS);
% S_Lrg_d_con = zeros(NX,DAYS);
% 
% S_Sml_f_nu = zeros(NX,DAYS);
% S_Sml_p_nu = zeros(NX,DAYS);
% S_Sml_d_nu = zeros(NX,DAYS);
% S_Med_f_nu = zeros(NX,DAYS);
% S_Med_p_nu = zeros(NX,DAYS);
% S_Med_d_nu = zeros(NX,DAYS);
% S_Lrg_p_nu = zeros(NX,DAYS);
% S_Lrg_d_nu = zeros(NX,DAYS);
% 
% S_Sml_f_prod = zeros(NX,DAYS);
% S_Sml_p_prod = zeros(NX,DAYS);
% S_Sml_d_prod = zeros(NX,DAYS);
% S_Med_f_prod = zeros(NX,DAYS);
% S_Med_p_prod = zeros(NX,DAYS);
% S_Med_d_prod = zeros(NX,DAYS);
% S_Lrg_p_prod = zeros(NX,DAYS);
% S_Lrg_d_prod = zeros(NX,DAYS);
% 
% S_Sml_f_gamma = zeros(NX,DAYS);
% S_Sml_p_gamma = zeros(NX,DAYS);
% S_Sml_d_gamma = zeros(NX,DAYS);
% S_Med_f_gamma = zeros(NX,DAYS);
% S_Med_p_gamma = zeros(NX,DAYS);
% S_Med_d_gamma = zeros(NX,DAYS);
% S_Lrg_p_gamma = zeros(NX,DAYS);
% S_Lrg_d_gamma = zeros(NX,DAYS);
% 
% S_Med_f_rep = zeros(NX,DAYS);
% S_Lrg_p_rep = zeros(NX,DAYS);
% S_Lrg_d_rep = zeros(NX,DAYS);
%
% S_Sml_f_die = zeros(NX,DAYS);
% S_Sml_p_die = zeros(NX,DAYS);
% S_Sml_d_die = zeros(NX,DAYS);
% S_Med_f_die = zeros(NX,DAYS);
% S_Med_p_die = zeros(NX,DAYS);
% S_Med_d_die = zeros(NX,DAYS);
% S_Lrg_p_die = zeros(NX,DAYS);
% S_Lrg_d_die = zeros(NX,DAYS);
% 
% S_Sml_f_clev = zeros(NX,DAYS);
% S_Sml_p_clev = zeros(NX,DAYS);
% S_Sml_d_clev = zeros(NX,DAYS);
% S_Med_f_clev = zeros(NX,DAYS);
% S_Med_p_clev = zeros(NX,DAYS);
% S_Med_d_clev = zeros(NX,DAYS);
% S_Lrg_p_clev = zeros(NX,DAYS);
% S_Lrg_d_clev = zeros(NX,DAYS);

%! Initialize
%!From a previous run
load([outdir 'Last_mo_Spinup_',vers,'_All_fishobs_',simname,'.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
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

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
% vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
% vidrecSF    = netcdf.defVar(ncidSF,'rec','double',[xy_dim,time_dim]);
% vidconSF    = netcdf.defVar(ncidSF,'con','double',[xy_dim,time_dim]);
% vidnuSF     = netcdf.defVar(ncidSF,'nu','double',[xy_dim,time_dim]);
% vidgammaSF  = netcdf.defVar(ncidSF,'gamma','double',[xy_dim,time_dim]);
% viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
% vidclevSF   = netcdf.defVar(ncidSF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
% vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
% vidrecSP    = netcdf.defVar(ncidSP,'rec','double',[xy_dim,time_dim]);
% vidconSP    = netcdf.defVar(ncidSP,'con','double',[xy_dim,time_dim]);
% vidnuSP     = netcdf.defVar(ncidSP,'nu','double',[xy_dim,time_dim]);
% vidgammaSP  = netcdf.defVar(ncidSP,'gamma','double',[xy_dim,time_dim]);
% viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
% vidclevSP   = netcdf.defVar(ncidSP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,time_dim]);
% vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
% vidrecSD    = netcdf.defVar(ncidSD,'rec','double',[xy_dim,time_dim]);
% vidconSD    = netcdf.defVar(ncidSD,'con','double',[xy_dim,time_dim]);
% vidnuSD     = netcdf.defVar(ncidSD,'nu','double',[xy_dim,time_dim]);
% vidgammaSD  = netcdf.defVar(ncidSD,'gamma','double',[xy_dim,time_dim]);
% viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
% vidclevSD   = netcdf.defVar(ncidSD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
vidfishMF   = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
% vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
% vidrecMF    = netcdf.defVar(ncidMF,'rec','double',[xy_dim,time_dim]);
% vidconMF    = netcdf.defVar(ncidMF,'con','double',[xy_dim,time_dim]);
% vidnuMF     = netcdf.defVar(ncidMF,'nu','double',[xy_dim,time_dim]);
% vidgammaMF  = netcdf.defVar(ncidMF,'gamma','double',[xy_dim,time_dim]);
% vidrepMF    = netcdf.defVar(ncidMF,'rep','double',[xy_dim,time_dim]);
% viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
% vidclevMF   = netcdf.defVar(ncidMF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
vidfishMP   = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
% vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
% vidrecMP    = netcdf.defVar(ncidMP,'rec','double',[xy_dim,time_dim]);
% vidconMP    = netcdf.defVar(ncidMP,'con','double',[xy_dim,time_dim]);
% vidnuMP     = netcdf.defVar(ncidMP,'nu','double',[xy_dim,time_dim]);
% vidgammaMP  = netcdf.defVar(ncidMP,'gamma','double',[xy_dim,time_dim]);
% viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
% vidclevMP   = netcdf.defVar(ncidMP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
vidfishMD   = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
% vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
% vidrecMD    = netcdf.defVar(ncidMD,'rec','double',[xy_dim,time_dim]);
% vidconMD    = netcdf.defVar(ncidMD,'con','double',[xy_dim,time_dim]);
% vidnuMD     = netcdf.defVar(ncidMD,'nu','double',[xy_dim,time_dim]);
% vidgammaMD  = netcdf.defVar(ncidMD,'gamma','double',[xy_dim,time_dim]);
% viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
% vidclevMD   = netcdf.defVar(ncidMD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
vidfishLP   = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
% vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
% vidrecLP    = netcdf.defVar(ncidLP,'rec','double',[xy_dim,time_dim]);
% vidconLP    = netcdf.defVar(ncidLP,'con','double',[xy_dim,time_dim]);
% vidnuLP     = netcdf.defVar(ncidLP,'nu','double',[xy_dim,time_dim]);
% vidgammaLP  = netcdf.defVar(ncidLP,'gamma','double',[xy_dim,time_dim]);
% vidrepLP    = netcdf.defVar(ncidLP,'rep','double',[xy_dim,time_dim]);
% viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
% vidclevLP   = netcdf.defVar(ncidLP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
vidfishLD   = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
% vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
% vidrecLD    = netcdf.defVar(ncidLD,'rec','double',[xy_dim,time_dim]);
% vidconLD    = netcdf.defVar(ncidLD,'con','double',[xy_dim,time_dim]);
% vidnuLD     = netcdf.defVar(ncidLD,'nu','double',[xy_dim,time_dim]);
% vidgammaLD  = netcdf.defVar(ncidLD,'gamma','double',[xy_dim,time_dim]);
% vidrepLD    = netcdf.defVar(ncidLD,'rep','double',[xy_dim,time_dim]);
% viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
% vidclevLD   = netcdf.defVar(ncidLD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

xy_dim     = netcdf.defDim(ncidB,'nid',NX);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);


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
        
        S_Bent_bio(:,DY) = BENT.mass;
        
        S_Sml_f(:,DY) = Sml_f.bio;
        S_Sml_p(:,DY) = Sml_p.bio;
        S_Sml_d(:,DY) = Sml_d.bio;
        S_Med_f(:,DY) = Med_f.bio;
        S_Med_p(:,DY) = Med_p.bio;
        S_Med_d(:,DY) = Med_d.bio;
        S_Lrg_p(:,DY) = Lrg_p.bio;
        S_Lrg_d(:,DY) = Lrg_d.bio;
        
        S_Med_f_fish(:,DY) = Med_f.caught;
        S_Med_p_fish(:,DY) = Med_p.caught;
        S_Med_d_fish(:,DY) = Med_d.caught;
        S_Lrg_p_fish(:,DY) = Lrg_p.caught;
        S_Lrg_d_fish(:,DY) = Lrg_d.caught;
        
        %         S_Sml_f_rec(:,DY) = Sml_f.rec;
        %         S_Sml_p_rec(:,DY) = Sml_p.rec;
        %         S_Sml_d_rec(:,DY) = Sml_d.rec;
        %         S_Med_f_rec(:,DY) = Med_f.rec;
        %         S_Med_p_rec(:,DY) = Med_p.rec;
        %         S_Med_d_rec(:,DY) = Med_d.rec;
        %         S_Lrg_p_rec(:,DY) = Lrg_p.rec;
        %         S_Lrg_d_rec(:,DY) = Lrg_d.rec;
        %
        %         S_Sml_f_con(:,DY) = Sml_f.I;
        %         S_Sml_p_con(:,DY) = Sml_p.I;
        %         S_Sml_d_con(:,DY) = Sml_d.I;
        %         S_Med_f_con(:,DY) = Med_f.I;
        %         S_Med_p_con(:,DY) = Med_p.I;
        %         S_Med_d_con(:,DY) = Med_d.I;
        %         S_Lrg_p_con(:,DY) = Lrg_p.I;
        %         S_Lrg_d_con(:,DY) = Lrg_d.I;
        %
        %         S_Sml_f_nu(:,DY) = Sml_f.nu;
        %         S_Sml_p_nu(:,DY) = Sml_p.nu;
        %         S_Sml_d_nu(:,DY) = Sml_d.nu;
        %         S_Med_f_nu(:,DY) = Med_f.nu;
        %         S_Med_p_nu(:,DY) = Med_p.nu;
        %         S_Med_d_nu(:,DY) = Med_d.nu;
        %         S_Lrg_p_nu(:,DY) = Lrg_p.nu;
        %         S_Lrg_d_nu(:,DY) = Lrg_d.nu;
        %
        %         S_Sml_f_prod(:,DY) = Sml_f.prod;
        %         S_Sml_p_prod(:,DY) = Sml_p.prod;
        %         S_Sml_d_prod(:,DY) = Sml_d.prod;
        %         S_Med_f_prod(:,DY) = Med_f.prod;
        %         S_Med_p_prod(:,DY) = Med_p.prod;
        %         S_Med_d_prod(:,DY) = Med_d.prod;
        %         S_Lrg_p_prod(:,DY) = Lrg_p.prod;
        %         S_Lrg_d_prod(:,DY) = Lrg_d.prod;
        %
        %         S_Sml_f_gamma(:,DY) = Sml_f.gamma;
        %         S_Sml_p_gamma(:,DY) = Sml_p.gamma;
        %         S_Sml_d_gamma(:,DY) = Sml_d.gamma;
        %         S_Med_f_gamma(:,DY) = Med_f.gamma;
        %         S_Med_p_gamma(:,DY) = Med_p.gamma;
        %         S_Med_d_gamma(:,DY) = Med_d.gamma;
        %         S_Lrg_p_gamma(:,DY) = Lrg_p.gamma;
        %         S_Lrg_d_gamma(:,DY) = Lrg_d.gamma;
        %
        %         S_Med_f_rep(:,DY) = Med_f.rep;
        %         S_Lrg_p_rep(:,DY) = Lrg_p.rep;
        %         S_Lrg_d_rep(:,DY) = Lrg_d.rep;
        %
        %         S_Sml_f_die(:,DY) = Sml_f.die;
        %         S_Sml_p_die(:,DY) = Sml_p.die;
        %         S_Sml_d_die(:,DY) = Sml_d.die;
        %         S_Med_f_die(:,DY) = Med_f.die;
        %         S_Med_p_die(:,DY) = Med_p.die;
        %         S_Med_d_die(:,DY) = Med_d.die;
        %         S_Lrg_p_die(:,DY) = Lrg_p.die;
        %         S_Lrg_d_die(:,DY) = Lrg_d.die;
        %
        %         S_Sml_f_clev(:,DY) = Sml_f.clev;
        %         S_Sml_p_clev(:,DY) = Sml_p.clev;
        %         S_Sml_d_clev(:,DY) = Sml_d.clev;
        %         S_Med_f_clev(:,DY) = Med_f.clev;
        %         S_Med_p_clev(:,DY) = Med_p.clev;
        %         S_Med_d_clev(:,DY) = Med_d.clev;
        %         S_Lrg_p_clev(:,DY) = Lrg_p.clev;
        %         S_Lrg_d_clev(:,DY) = Lrg_d.clev;
        
    end %Days
    
    %! Calculate monthly means and save
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
        
        netcdf.putVar(ncidMF,vidfishMF,[0 MNT-1],[NX 1],mean(S_Med_f_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidfishMP,[0 MNT-1],[NX 1],mean(S_Med_p_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidfishMD,[0 MNT-1],[NX 1],mean(S_Med_d_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidfishLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidfishLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_fish(:,a(i):b(i)),2));
        
        %           netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(S_Sml_f_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(S_Sml_p_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,vidprodSD,[0 MNT-1],[NX 1],mean(S_Sml_d_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(S_Med_f_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(S_Med_p_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,vidprodMD,[0 MNT-1],[NX 1],mean(S_Med_d_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_prod(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidprodLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_prod(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidSF,vidrecSF,[0 MNT-1],[NX 1],mean(S_Sml_f_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,vidrecSP,[0 MNT-1],[NX 1],mean(S_Sml_p_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,vidrecSD,[0 MNT-1],[NX 1],mean(S_Sml_d_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,vidrecMF,[0 MNT-1],[NX 1],mean(S_Med_f_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,vidrecMP,[0 MNT-1],[NX 1],mean(S_Med_p_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,vidrecMD,[0 MNT-1],[NX 1],mean(S_Med_d_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidrecLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rec(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidrecLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rec(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidSF,vidconSF,[0 MNT-1],[NX 1],mean(S_Sml_f_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,vidconSP,[0 MNT-1],[NX 1],mean(S_Sml_p_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,vidconSD,[0 MNT-1],[NX 1],mean(S_Sml_d_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,vidconMF,[0 MNT-1],[NX 1],mean(S_Med_f_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,vidconMP,[0 MNT-1],[NX 1],mean(S_Med_p_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,vidconMD,[0 MNT-1],[NX 1],mean(S_Med_d_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidconLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_con(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidconLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_con(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidSF,vidnuSF,[0 MNT-1],[NX 1],mean(S_Sml_f_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,vidnuSP,[0 MNT-1],[NX 1],mean(S_Sml_p_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,vidnuSD,[0 MNT-1],[NX 1],mean(S_Sml_d_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,vidnuMF,[0 MNT-1],[NX 1],mean(S_Med_f_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,vidnuMP,[0 MNT-1],[NX 1],mean(S_Med_p_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,vidnuMD,[0 MNT-1],[NX 1],mean(S_Med_d_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidnuLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_nu(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidnuLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_nu(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidSF,vidgammaSF,[0 MNT-1],[NX 1],mean(S_Sml_f_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,vidgammaSP,[0 MNT-1],[NX 1],mean(S_Sml_p_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,vidgammaSD,[0 MNT-1],[NX 1],mean(S_Sml_d_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,vidgammaMF,[0 MNT-1],[NX 1],mean(S_Med_f_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,vidgammaMP,[0 MNT-1],[NX 1],mean(S_Med_p_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,vidgammaMD,[0 MNT-1],[NX 1],mean(S_Med_d_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidgammaLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_gamma(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidgammaLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_gamma(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidSF,vidclevSF,[0 MNT-1],[NX 1],mean(S_Sml_f_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,vidclevSP,[0 MNT-1],[NX 1],mean(S_Sml_p_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,vidclevSD,[0 MNT-1],[NX 1],mean(S_Sml_d_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,vidclevMF,[0 MNT-1],[NX 1],mean(S_Med_f_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,vidclevMP,[0 MNT-1],[NX 1],mean(S_Med_p_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,vidclevMD,[0 MNT-1],[NX 1],mean(S_Med_d_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidclevLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_clev(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidclevLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_clev(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidSF,viddieSF,[0 MNT-1],[NX 1],mean(S_Sml_f_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSP,viddieSP,[0 MNT-1],[NX 1],mean(S_Sml_p_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidSD,viddieSD,[0 MNT-1],[NX 1],mean(S_Sml_d_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMF,viddieMF,[0 MNT-1],[NX 1],mean(S_Med_f_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMP,viddieMP,[0 MNT-1],[NX 1],mean(S_Med_p_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidMD,viddieMD,[0 MNT-1],[NX 1],mean(S_Med_d_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,viddieLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_die(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,viddieLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_die(:,a(i):b(i)),2));
        %
        %           netcdf.putVar(ncidMF,vidrepMF,[0 MNT-1],[NX 1],mean(S_Med_f_rep(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLP,vidrepLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rep(:,a(i):b(i)),2));
        % 			netcdf.putVar(ncidLD,vidrepLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rep(:,a(i):b(i)),2));
        %
        
    end %Monthly mean
    
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

end %esms

end
