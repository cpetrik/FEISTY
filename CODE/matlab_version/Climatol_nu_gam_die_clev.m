%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_nu_gam_die_clev()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit h gam kt bpow
global bent_eff rfrac CC D J Sm A benc bcmx
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn
global tstep K CGRD ni nj

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
frate = 0.3; %Fish(F);
dfrate = frate/365.0;
%0=no fishing; 1=fishing
if (frate>0)
    harv = 1;
else
    harv = 0;
end


%! Make core parameters/constants (global)
make_parameters()

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
Pdrpbx = '/Users/cpetrik/Dropbox/';
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat');
NX = length(GRD.Z);
ID = 1:NX;
%2D Grid for advect-diff
load([Pdrpbx 'Princeton/POEM_2.0/CODE/Data/Hindcast_cgrid_cp2D.mat']);
[ni,nj] = size(CGRD.mask);

%! Create a directory for output
tfcrit = num2str(int64(100*fcrit));
tkap = num2str(100+int64(10*K_a));
td = num2str(1000+int64(100*LD_phi_MP));
tj = num2str(1000+int64(100*MP_phi_S));
tsm = num2str(1000+int64(100*MF_phi_MZ));
ta = num2str(1000+int64(100*LP_phi_MF));
tbe = num2str(100+int64(100*bent_eff));
tmort = num2str(MORT);
tcc = num2str(1000+int64(100*CC));
tre = num2str(100000+int64(round(10000*rfrac)));
tre2 = num2str(100000+int64(round(10000*rfrac*1)));
if (frate >= 0.1)
    tfish = num2str(100+int64(10*frate));
    tF = num2str(1000+int64(100*frate*MFsel));
    tP = num2str(1000+int64(100*frate*LPsel));
    tD = num2str(1000+int64(100*frate*LDsel));
    tJ = num2str(100+int64(10*Jsel));
else
    tfish = num2str(1000+int64(100*frate));
    tF = num2str(1000+int64(100*frate*MFsel));
    tP = num2str(1000+int64(100*frate*LPsel));
    tD = num2str(1000+int64(100*frate*LDsel));
end
if (MFsel > 0)
    if (LPsel > 0 && LDsel > 0)
        sel='All';
    else
        sel='F';
    end
else
    if (LPsel > 0 && LDsel > 0)
        sel = 'L';
    elseif (LPsel > 0)
        sel = 'P';
    elseif (LDsel > 0)
        sel = 'D';
    end
end
if (pdc == 0)
    coup = 'NoDc';
elseif (pdc == 1)
    coup = 'Dc';
elseif (pdc == 2)
    coup = 'PDc';
end
tcfn = num2str(h);
tefn = num2str(round(gam));
tkfn = num2str(100+int64(100*kt));
tbfn = num2str(1000+int64(1000*bpow));
tbenc = num2str(1000+int64(1000*benc));
tbcmx = num2str(1000+int64(1000*bcmx));
%simname = [coup,'_enc',tefn,'_cmax-metab',tcfn,'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
%simname = [coup,'_enc',tefn,'_cmax-metab',tcfn,'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
%simname = ['Diff_',coup,'_enc',tefn,'_cmax-metab',tcfn,'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
%simname = [coup,'_enc',tefn,'_cmax-metab',tcfn,'_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
%simname = [coup,'_enc',tefn,'_cm',tcfn,'_m-b200_c-b',tbfn(2:end),'_m-k',tkfn(2:end),'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
%simname = [coup,'_enc',tefn,'-b',tbfn(2:end),'_cm',tcfn,'_m-b175-k',tkfn(2:end),'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
%simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_cm',tcfn,'_m-b',tbfn(2:end),'-k',tkfn(2:end),'_fcrit',tfcrit,'_c-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_cm',tcfn,'_m-b',tbfn(2:end),'-k',tkfn(2:end),'_fcrit',tfcrit,'_c-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
if (~isdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname]))
    mkdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname])
end
% if (~isdir([Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/',simname]))
%     mkdir([Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/',simname])
% end

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

S_Sml_f_clev = zeros(NX,DAYS);
S_Sml_p_clev = zeros(NX,DAYS);
S_Sml_d_clev = zeros(NX,DAYS);
S_Med_f_clev = zeros(NX,DAYS);
S_Med_p_clev = zeros(NX,DAYS);
S_Med_d_clev = zeros(NX,DAYS);
S_Lrg_p_clev = zeros(NX,DAYS);
S_Lrg_d_clev = zeros(NX,DAYS);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
Med_d.td(1:NX) = 0.0;
Lrg_d.td(1:NX) = 0.0;
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
if harv==0
    file_sml_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_sml_f.nc'];
    file_sml_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_sml_p.nc'];
    file_sml_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_sml_d.nc'];
    file_med_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_med_f.nc'];
    file_med_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_med_p.nc'];
    file_med_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_med_d.nc'];
    file_lrg_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_lrg_p.nc'];
    file_lrg_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_nu_gam_die_clev_lrg_d.nc'];
else
    file_sml_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_sml_f.nc'];
    file_sml_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_sml_p.nc'];
    file_sml_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_sml_d.nc'];
    file_med_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_med_f.nc'];
    file_med_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_med_p.nc'];
    file_med_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_med_d.nc'];
    file_lrg_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_lrg_p.nc'];
    file_lrg_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_nu_gam_die_clev_lrg_d.nc'];
    
%     file_sml_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_sml_f.nc'];
%     file_sml_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_sml_p.nc'];
%     file_sml_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_sml_d.nc'];
%     file_med_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_med_f.nc'];
%     file_med_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_med_p.nc'];
%     file_med_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_med_d.nc'];
%     file_lrg_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_lrg_p.nc'];
%     file_lrg_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_Juve',tJ(2:end),'_lrg_d.nc'];
end

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
%oldFormat = netcdf.setDefaultFormat('NC_FORMAT_64BIT');
netcdf.setDefaultFormat('NC_FORMAT_64BIT')

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidnuSF     = netcdf.defVar(ncidSF,'nu','double',[xy_dim,time_dim]);
vidgammaSF  = netcdf.defVar(ncidSF,'gamma','double',[xy_dim,time_dim]);
viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
vidclevSF   = netcdf.defVar(ncidSF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidnuSP     = netcdf.defVar(ncidSP,'nu','double',[xy_dim,time_dim]);
vidgammaSP  = netcdf.defVar(ncidSP,'gamma','double',[xy_dim,time_dim]);
viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
vidclevSP   = netcdf.defVar(ncidSP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidnuSD     = netcdf.defVar(ncidSD,'nu','double',[xy_dim,time_dim]);
vidgammaSD  = netcdf.defVar(ncidSD,'gamma','double',[xy_dim,time_dim]);
viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
vidclevSD   = netcdf.defVar(ncidSD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidnuMF     = netcdf.defVar(ncidMF,'nu','double',[xy_dim,time_dim]);
vidgammaMF  = netcdf.defVar(ncidMF,'gamma','double',[xy_dim,time_dim]);
viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
vidclevMF   = netcdf.defVar(ncidMF,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidnuMP     = netcdf.defVar(ncidMP,'nu','double',[xy_dim,time_dim]);
vidgammaMP  = netcdf.defVar(ncidMP,'gamma','double',[xy_dim,time_dim]);
viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
vidclevMP   = netcdf.defVar(ncidMP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidnuMD     = netcdf.defVar(ncidMD,'nu','double',[xy_dim,time_dim]);
vidgammaMD  = netcdf.defVar(ncidMD,'gamma','double',[xy_dim,time_dim]);
viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
vidclevMD   = netcdf.defVar(ncidMD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidnuLP     = netcdf.defVar(ncidLP,'nu','double',[xy_dim,time_dim]);
vidgammaLP  = netcdf.defVar(ncidLP,'gamma','double',[xy_dim,time_dim]);
viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
vidclevLP   = netcdf.defVar(ncidLP,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidnuLD     = netcdf.defVar(ncidLD,'nu','double',[xy_dim,time_dim]);
vidgammaLD  = netcdf.defVar(ncidLD,'gamma','double',[xy_dim,time_dim]);
viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
vidclevLD   = netcdf.defVar(ncidLD,'clev','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Run model with no fishing
MNT=0;
for YR = 1:YEARS % years
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate,CC);
        
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
        
        S_Sml_f_clev(:,DY) = Sml_f.clev;
        S_Sml_p_clev(:,DY) = Sml_p.clev;
        S_Sml_d_clev(:,DY) = Sml_d.clev;
        S_Med_f_clev(:,DY) = Med_f.clev;
        S_Med_p_clev(:,DY) = Med_p.clev;
        S_Med_d_clev(:,DY) = Med_d.clev;
        S_Lrg_p_clev(:,DY) = Lrg_p.clev;
        S_Lrg_d_clev(:,DY) = Lrg_d.clev;
        
        
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
        
        netcdf.putVar(ncidSF,vidclevSF,[0 MNT-1],[NX 1],mean(S_Sml_f_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidclevSP,[0 MNT-1],[NX 1],mean(S_Sml_p_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidclevSD,[0 MNT-1],[NX 1],mean(S_Sml_d_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidclevMF,[0 MNT-1],[NX 1],mean(S_Med_f_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidclevMP,[0 MNT-1],[NX 1],mean(S_Med_p_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidclevMD,[0 MNT-1],[NX 1],mean(S_Med_d_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidclevLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_clev(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidclevLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_clev(:,a(i):b(i)),2));
        
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