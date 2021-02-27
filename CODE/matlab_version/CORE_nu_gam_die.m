%%%%!! RUN CORE-FORCED HISTORIC WITH FISHING FOR ALL LOCATIONS
function CORE_nu_gam_die()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
frate = 0.3;
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_parameters() % make core parameters/constants

%! Grid
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

%! How long to run the model
YEARS = length(1950:2007);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
[fname,simname] = sub_fname_core(frate);

%! Storage variables
S_Sml_f_nu = zeros(NX,DAYS);
S_Sml_p_nu = zeros(NX,DAYS);
S_Sml_d_nu = zeros(NX,DAYS);
S_Med_f_nu = zeros(NX,DAYS);
S_Med_p_nu = zeros(NX,DAYS);
S_Med_d_nu = zeros(NX,DAYS);
S_Lrg_p_nu = zeros(NX,DAYS);
S_Lrg_d_nu = zeros(NX,DAYS);

S_Sml_f_gamma = zeros(NX,DAYS);
S_Sml_p_gamma = zeros(NX,DAYS);
S_Sml_d_gamma = zeros(NX,DAYS);
S_Med_f_gamma = zeros(NX,DAYS);
S_Med_p_gamma = zeros(NX,DAYS);
S_Med_d_gamma = zeros(NX,DAYS);
S_Lrg_p_gamma = zeros(NX,DAYS);
S_Lrg_d_gamma = zeros(NX,DAYS);

S_Sml_f_die = zeros(NX,DAYS);
S_Sml_p_die = zeros(NX,DAYS);
S_Sml_d_die = zeros(NX,DAYS);
S_Med_f_die = zeros(NX,DAYS);
S_Med_p_die = zeros(NX,DAYS);
S_Med_d_die = zeros(NX,DAYS);
S_Lrg_p_die = zeros(NX,DAYS);
S_Lrg_d_die = zeros(NX,DAYS);

%% ! Initialize
init_sim = simname;
load(['/Volumes/FEISTY/NC/Matlab_new_size/',init_sim '/Last_mo_spinup_fished_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_nu_gam_die_sml_f.nc'];
file_sml_p = [fname,'_nu_gam_die_sml_p.nc'];
file_sml_d = [fname,'_nu_gam_die_sml_d.nc'];
file_med_f = [fname,'_nu_gam_die_med_f.nc'];
file_med_p = [fname,'_nu_gam_die_med_p.nc'];
file_med_d = [fname,'_nu_gam_die_med_d.nc'];
file_lrg_p = [fname,'_nu_gam_die_lrg_p.nc'];
file_lrg_d = [fname,'_nu_gam_die_lrg_d.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidnuSF     = netcdf.defVar(ncidSF,'nu','double',[xy_dim,time_dim]);
vidgammaSF  = netcdf.defVar(ncidSF,'gamma','double',[xy_dim,time_dim]);
viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidnuSP     = netcdf.defVar(ncidSP,'nu','double',[xy_dim,time_dim]);
vidgammaSP  = netcdf.defVar(ncidSP,'gamma','double',[xy_dim,time_dim]);
viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidnuSD     = netcdf.defVar(ncidSD,'nu','double',[xy_dim,time_dim]);
vidgammaSD  = netcdf.defVar(ncidSD,'gamma','double',[xy_dim,time_dim]);
viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidnuMF     = netcdf.defVar(ncidMF,'nu','double',[xy_dim,time_dim]);
vidgammaMF  = netcdf.defVar(ncidMF,'gamma','double',[xy_dim,time_dim]);
viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidnuMP     = netcdf.defVar(ncidMP,'nu','double',[xy_dim,time_dim]);
vidgammaMP  = netcdf.defVar(ncidMP,'gamma','double',[xy_dim,time_dim]);
viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidnuMD     = netcdf.defVar(ncidMD,'nu','double',[xy_dim,time_dim]);
vidgammaMD  = netcdf.defVar(ncidMD,'gamma','double',[xy_dim,time_dim]);
viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidnuLP     = netcdf.defVar(ncidLP,'nu','double',[xy_dim,time_dim]);
vidgammaLP  = netcdf.defVar(ncidLP,'gamma','double',[xy_dim,time_dim]);
viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidnuLD     = netcdf.defVar(ncidLD,'nu','double',[xy_dim,time_dim]);
vidgammaLD  = netcdf.defVar(ncidLD,'gamma','double',[xy_dim,time_dim]);
viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    %! Load a year's COBALT data
    ti = num2str(YR+1949);
    load(['/Volumes/MIP/GCM_DATA/CORE-forced/Data_ocean_cobalt_daily_',ti,'.mat'],'COBALT');
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ',num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
        
        %! Store
        S_Sml_f_nu(:,DY) = Sml_f.nu;
        S_Sml_p_nu(:,DY) = Sml_p.nu;
        S_Sml_d_nu(:,DY) = Sml_d.nu;
        S_Med_f_nu(:,DY) = Med_f.nu;
        S_Med_p_nu(:,DY) = Med_p.nu;
        S_Med_d_nu(:,DY) = Med_d.nu;
        S_Lrg_p_nu(:,DY) = Lrg_p.nu;
        S_Lrg_d_nu(:,DY) = Lrg_d.nu;
        
        S_Sml_f_gamma(:,DY) = Sml_f.gamma;
        S_Sml_p_gamma(:,DY) = Sml_p.gamma;
        S_Sml_d_gamma(:,DY) = Sml_d.gamma;
        S_Med_f_gamma(:,DY) = Med_f.gamma;
        S_Med_p_gamma(:,DY) = Med_p.gamma;
        S_Med_d_gamma(:,DY) = Med_d.gamma;
        S_Lrg_p_gamma(:,DY) = Lrg_p.gamma;
        S_Lrg_d_gamma(:,DY) = Lrg_d.gamma;
               
        S_Sml_f_die(:,DY) = Sml_f.die;
        S_Sml_p_die(:,DY) = Sml_p.die;
        S_Sml_d_die(:,DY) = Sml_d.die;
        S_Med_f_die(:,DY) = Med_f.die;
        S_Med_p_die(:,DY) = Med_p.die;
        S_Med_d_die(:,DY) = Med_d.die;
        S_Lrg_p_die(:,DY) = Lrg_p.die;
        S_Lrg_d_die(:,DY) = Lrg_d.die;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidnuSF,[0 MNT-1],[NX 1],mean(S_Sml_f_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidnuSP,[0 MNT-1],[NX 1],mean(S_Sml_p_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidnuSD,[0 MNT-1],[NX 1],mean(S_Sml_d_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidnuMF,[0 MNT-1],[NX 1],mean(S_Med_f_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidnuMP,[0 MNT-1],[NX 1],mean(S_Med_p_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidnuMD,[0 MNT-1],[NX 1],mean(S_Med_d_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidnuLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_nu(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidnuLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_nu(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidSF,vidgammaSF,[0 MNT-1],[NX 1],mean(S_Sml_f_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidgammaSP,[0 MNT-1],[NX 1],mean(S_Sml_p_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidgammaSD,[0 MNT-1],[NX 1],mean(S_Sml_d_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidgammaMF,[0 MNT-1],[NX 1],mean(S_Med_f_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidgammaMP,[0 MNT-1],[NX 1],mean(S_Med_p_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidgammaMD,[0 MNT-1],[NX 1],mean(S_Med_d_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidgammaLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_gamma(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidgammaLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_gamma(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidSF,viddieSF,[0 MNT-1],[NX 1],mean(S_Sml_f_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,viddieSP,[0 MNT-1],[NX 1],mean(S_Sml_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,viddieSD,[0 MNT-1],[NX 1],mean(S_Sml_d_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,viddieMF,[0 MNT-1],[NX 1],mean(S_Med_f_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,viddieMP,[0 MNT-1],[NX 1],mean(S_Med_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,viddieMD,[0 MNT-1],[NX 1],mean(S_Med_d_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,viddieLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,viddieLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_die(:,a(i):b(i)),2));
              
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


end
