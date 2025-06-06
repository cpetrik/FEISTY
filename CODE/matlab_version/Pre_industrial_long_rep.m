%%%%!! RUN PRE-INDUSTRIAL FOR ALL LOCATIONS
function Pre_industrial_long_rep()

% 400-yr preindustrial

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
frate = 0;
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_parameters() % make core parameters/constants

%! Grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

%! How long to run the model
YEARS = 400;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
fname = sub_fname_pi400(frate);

%! Storage variables
% S_Bent_bio = zeros(NX,DAYS);
% 
% S_Sml_f = zeros(NX,DAYS);
% S_Sml_p = zeros(NX,DAYS);
% S_Sml_d = zeros(NX,DAYS);
% S_Med_f = zeros(NX,DAYS);
% S_Med_p = zeros(NX,DAYS);
% S_Med_d = zeros(NX,DAYS);
% S_Lrg_p = zeros(NX,DAYS);
% S_Lrg_d = zeros(NX,DAYS);

% S_Sml_f_rec = zeros(NX,DAYS);
% S_Sml_p_rec = zeros(NX,DAYS);
% S_Sml_d_rec = zeros(NX,DAYS);
% S_Med_f_rec = zeros(NX,DAYS);
% S_Med_p_rec = zeros(NX,DAYS);
% S_Med_d_rec = zeros(NX,DAYS);
% S_Lrg_p_rec = zeros(NX,DAYS);
% S_Lrg_d_rec = zeros(NX,DAYS);

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

% S_Sml_f_prod = zeros(NX,DAYS);
% S_Sml_p_prod = zeros(NX,DAYS);
% S_Sml_d_prod = zeros(NX,DAYS);
% S_Med_f_prod = zeros(NX,DAYS);
% S_Med_p_prod = zeros(NX,DAYS);
% S_Med_d_prod = zeros(NX,DAYS);
% S_Lrg_p_prod = zeros(NX,DAYS);
% S_Lrg_d_prod = zeros(NX,DAYS);

% S_Sml_f_gamma = zeros(NX,DAYS);
% S_Sml_p_gamma = zeros(NX,DAYS);
% S_Sml_d_gamma = zeros(NX,DAYS);
% S_Med_f_gamma = zeros(NX,DAYS);
% S_Med_p_gamma = zeros(NX,DAYS);
% S_Med_d_gamma = zeros(NX,DAYS);
% S_Lrg_p_gamma = zeros(NX,DAYS);
% S_Lrg_d_gamma = zeros(NX,DAYS);

S_Med_f_rep = zeros(NX,DAYS);
S_Lrg_p_rep = zeros(NX,DAYS);
S_Lrg_d_rep = zeros(NX,DAYS);

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
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_med_f = [fname,'_rep_med_f.nc'];
file_lrg_p = [fname,'_rep_lrg_p.nc'];
file_lrg_d = [fname,'_rep_lrg_d.nc'];

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
% vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
% vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
% vidrecMF    = netcdf.defVar(ncidMF,'rec','double',[xy_dim,time_dim]);
% vidconMF    = netcdf.defVar(ncidMF,'con','double',[xy_dim,time_dim]);
% vidnuMF     = netcdf.defVar(ncidMF,'nu','double',[xy_dim,time_dim]);
% vidgammaMF  = netcdf.defVar(ncidMF,'gamma','double',[xy_dim,time_dim]);
vidrepMF    = netcdf.defVar(ncidMF,'rep','double',[xy_dim,time_dim]);
% viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
% vidclevMF   = netcdf.defVar(ncidMF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
% vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
% vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
% vidrecLP    = netcdf.defVar(ncidLP,'rec','double',[xy_dim,time_dim]);
% vidconLP    = netcdf.defVar(ncidLP,'con','double',[xy_dim,time_dim]);
% vidnuLP     = netcdf.defVar(ncidLP,'nu','double',[xy_dim,time_dim]);
% vidgammaLP  = netcdf.defVar(ncidLP,'gamma','double',[xy_dim,time_dim]);
vidrepLP    = netcdf.defVar(ncidLP,'rep','double',[xy_dim,time_dim]);
% viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
% vidclevLP   = netcdf.defVar(ncidLP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
% vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
% vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
% vidrecLD    = netcdf.defVar(ncidLD,'rec','double',[xy_dim,time_dim]);
% vidconLD    = netcdf.defVar(ncidLD,'con','double',[xy_dim,time_dim]);
% vidnuLD     = netcdf.defVar(ncidLD,'nu','double',[xy_dim,time_dim]);
% vidgammaLD  = netcdf.defVar(ncidLD,'gamma','double',[xy_dim,time_dim]);
vidrepLD    = netcdf.defVar(ncidLD,'rep','double',[xy_dim,time_dim]);
% viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
% vidclevLD   = netcdf.defVar(ncidLD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    num2str(YR)
    %! Load a year's COBALT data
    ti = num2str(YR+1000);
    load(['/Volumes/MIP/GCM_DATA/ESM2M_PI/Data_ocean_cobalt_daily_',ti,'.mat']);
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
        
        %! Store
%         S_Bent_bio(:,DY) = BENT.mass;
%         
%         S_Sml_f(:,DY) = Sml_f.bio;
%         S_Sml_p(:,DY) = Sml_p.bio;
%         S_Sml_d(:,DY) = Sml_d.bio;
%         S_Med_f(:,DY) = Med_f.bio;
%         S_Med_p(:,DY) = Med_p.bio;
%         S_Med_d(:,DY) = Med_d.bio;
%         S_Lrg_p(:,DY) = Lrg_p.bio;
%         S_Lrg_d(:,DY) = Lrg_d.bio;
        
%         S_Sml_f_rec(:,DY) = Sml_f.rec;
%         S_Sml_p_rec(:,DY) = Sml_p.rec;
%         S_Sml_d_rec(:,DY) = Sml_d.rec;
%         S_Med_f_rec(:,DY) = Med_f.rec;
%         S_Med_p_rec(:,DY) = Med_p.rec;
%         S_Med_d_rec(:,DY) = Med_d.rec;
%         S_Lrg_p_rec(:,DY) = Lrg_p.rec;
%         S_Lrg_d_rec(:,DY) = Lrg_d.rec;
        
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
        
%         S_Sml_f_prod(:,DY) = Sml_f.prod;
%         S_Sml_p_prod(:,DY) = Sml_p.prod;
%         S_Sml_d_prod(:,DY) = Sml_d.prod;
%         S_Med_f_prod(:,DY) = Med_f.prod;
%         S_Med_p_prod(:,DY) = Med_p.prod;
%         S_Med_d_prod(:,DY) = Med_d.prod;
%         S_Lrg_p_prod(:,DY) = Lrg_p.prod;
%         S_Lrg_d_prod(:,DY) = Lrg_d.prod;
        
%         S_Sml_f_gamma(:,DY) = Sml_f.gamma;
%         S_Sml_p_gamma(:,DY) = Sml_p.gamma;
%         S_Sml_d_gamma(:,DY) = Sml_d.gamma;
%         S_Med_f_gamma(:,DY) = Med_f.gamma;
%         S_Med_p_gamma(:,DY) = Med_p.gamma;
%         S_Med_d_gamma(:,DY) = Med_d.gamma;
%         S_Lrg_p_gamma(:,DY) = Lrg_p.gamma;
%         S_Lrg_d_gamma(:,DY) = Lrg_d.gamma;
%         
        S_Med_f_rep(:,DY) = Med_f.rep;
        S_Lrg_p_rep(:,DY) = Lrg_p.rep;
        S_Lrg_d_rep(:,DY) = Lrg_d.rep;
        
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
%         netcdf.putVar(ncidB,vidbioB,[0 MNT-1],[NX 1],mean(S_Bent_bio(:,a(i):b(i)),2));
%         netcdf.putVar(ncidB,vidTB,MNT-1,1,MNT);
%         
%         netcdf.putVar(ncidSF,vidbioSF,[0 MNT-1],[NX 1],mean(S_Sml_f(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidbioSP,[0 MNT-1],[NX 1],mean(S_Sml_p(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidbioSD,[0 MNT-1],[NX 1],mean(S_Sml_d(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidbioMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidbioMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidbioLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidbioLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));
%         
%         netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(S_Sml_f_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(S_Sml_p_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidprodSD,[0 MNT-1],[NX 1],mean(S_Sml_d_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(S_Med_f_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(S_Med_p_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidprodMD,[0 MNT-1],[NX 1],mean(S_Med_d_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidprodLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_prod(:,a(i):b(i)),2));
%         
%         netcdf.putVar(ncidSF,vidrecSF,[0 MNT-1],[NX 1],mean(S_Sml_f_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidrecSP,[0 MNT-1],[NX 1],mean(S_Sml_p_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidrecSD,[0 MNT-1],[NX 1],mean(S_Sml_d_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidrecMF,[0 MNT-1],[NX 1],mean(S_Med_f_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidrecMP,[0 MNT-1],[NX 1],mean(S_Med_p_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidrecMD,[0 MNT-1],[NX 1],mean(S_Med_d_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidrecLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rec(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidrecLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rec(:,a(i):b(i)),2));
        
%         netcdf.putVar(ncidSF,vidconSF,[0 MNT-1],[NX 1],mean(S_Sml_f_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidconSP,[0 MNT-1],[NX 1],mean(S_Sml_p_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidconSD,[0 MNT-1],[NX 1],mean(S_Sml_d_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidconMF,[0 MNT-1],[NX 1],mean(S_Med_f_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidconMP,[0 MNT-1],[NX 1],mean(S_Med_p_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidconMD,[0 MNT-1],[NX 1],mean(S_Med_d_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidconLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_con(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidconLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_con(:,a(i):b(i)),2));
%         
%         netcdf.putVar(ncidSF,vidnuSF,[0 MNT-1],[NX 1],mean(S_Sml_f_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidnuSP,[0 MNT-1],[NX 1],mean(S_Sml_p_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidnuSD,[0 MNT-1],[NX 1],mean(S_Sml_d_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidnuMF,[0 MNT-1],[NX 1],mean(S_Med_f_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidnuMP,[0 MNT-1],[NX 1],mean(S_Med_p_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidnuMD,[0 MNT-1],[NX 1],mean(S_Med_d_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidnuLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_nu(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidnuLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_nu(:,a(i):b(i)),2));
%         
%         netcdf.putVar(ncidSF,vidgammaSF,[0 MNT-1],[NX 1],mean(S_Sml_f_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidgammaSP,[0 MNT-1],[NX 1],mean(S_Sml_p_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidgammaSD,[0 MNT-1],[NX 1],mean(S_Sml_d_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidgammaMF,[0 MNT-1],[NX 1],mean(S_Med_f_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidgammaMP,[0 MNT-1],[NX 1],mean(S_Med_p_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidgammaMD,[0 MNT-1],[NX 1],mean(S_Med_d_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidgammaLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_gamma(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidgammaLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_gamma(:,a(i):b(i)),2));
%         
%         netcdf.putVar(ncidSF,vidclevSF,[0 MNT-1],[NX 1],mean(S_Sml_f_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidclevSP,[0 MNT-1],[NX 1],mean(S_Sml_p_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,vidclevSD,[0 MNT-1],[NX 1],mean(S_Sml_d_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidclevMF,[0 MNT-1],[NX 1],mean(S_Med_f_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidclevMP,[0 MNT-1],[NX 1],mean(S_Med_p_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,vidclevMD,[0 MNT-1],[NX 1],mean(S_Med_d_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidclevLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_clev(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,vidclevLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_clev(:,a(i):b(i)),2));
%         
%         netcdf.putVar(ncidSF,viddieSF,[0 MNT-1],[NX 1],mean(S_Sml_f_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,viddieSP,[0 MNT-1],[NX 1],mean(S_Sml_p_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSD,viddieSD,[0 MNT-1],[NX 1],mean(S_Sml_d_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,viddieMF,[0 MNT-1],[NX 1],mean(S_Med_f_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,viddieMP,[0 MNT-1],[NX 1],mean(S_Med_p_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMD,viddieMD,[0 MNT-1],[NX 1],mean(S_Med_d_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,viddieLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_die(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLD,viddieLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_die(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidMF,vidrepMF,[0 MNT-1],[NX 1],mean(S_Med_f_rep(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidrepLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rep(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidrepLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rep(:,a(i):b(i)),2));
         
    end %Monthly mean
    
end %Years


%! Close save
netcdf.close(ncidMF);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
