%%%%!! RUN Climatol FOR ALL LOCATIONS
%%% Need consumption of each type by each type
function Climatol_con_types()

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
simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_cm',tcfn,'_m-b',tbfn(2:end),'-k',tkfn(2:end),'_fcrit',tfcrit,'_c-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_noCC_RE',tre(2:end)];
if (~isdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname]))
    mkdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname])
end

%! Storage variables
S_Sml_f_con = zeros(NX,DAYS);
S_Sml_p_con = zeros(NX,DAYS);
S_Sml_d_con = zeros(NX,DAYS);
S_Med_f_conZm = zeros(NX,DAYS);
S_Med_f_conZl = zeros(NX,DAYS);
S_Med_f_conF = zeros(NX,DAYS);
S_Med_f_conP = zeros(NX,DAYS);
S_Med_p_conZm = zeros(NX,DAYS);
S_Med_p_conZl = zeros(NX,DAYS);
S_Med_p_conF = zeros(NX,DAYS);
S_Med_p_conP = zeros(NX,DAYS);
S_Med_d_con = zeros(NX,DAYS);
S_Lrg_p_conF = zeros(NX,DAYS);
S_Lrg_p_conP = zeros(NX,DAYS);
S_Lrg_d_conB = zeros(NX,DAYS);
S_Lrg_d_conF = zeros(NX,DAYS);
S_Lrg_d_conP = zeros(NX,DAYS);
S_Lrg_d_conD = zeros(NX,DAYS);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
Med_d.td(1:NX) = 0.0;
Lrg_d.td(1:NX) = 0.0;
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
if harv==0
    file_sml_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_sml_f.nc'];
    file_sml_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_sml_p.nc'];
    file_sml_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_sml_d.nc'];
    file_med_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_med_f.nc'];
    file_med_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_med_p.nc'];
    file_med_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_med_d.nc'];
    file_lrg_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_lrg_p.nc'];
    file_lrg_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_pristine_con_types_lrg_d.nc'];
else
    file_sml_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_sml_f.nc'];
    file_sml_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_sml_p.nc'];
    file_sml_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_sml_d.nc'];
    file_med_f = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_med_f.nc'];
    file_med_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_med_p.nc'];
    file_med_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_med_d.nc'];
    file_lrg_p = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_lrg_p.nc'];
    file_lrg_d = ['/Volumes/GFDL/NC/Matlab_new_size/',simname, '/Climatol_', sel,'_fish',tfish(2:end),'_con_types_lrg_d.nc'];
    
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
vidconSF    = netcdf.defVar(ncidSF,'conZ','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidconSP    = netcdf.defVar(ncidSP,'conZ','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidconSD    = netcdf.defVar(ncidSD,'conZ','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidconMFzm  = netcdf.defVar(ncidMF,'conZm','double',[xy_dim,time_dim]);
vidconMFzl  = netcdf.defVar(ncidMF,'conZl','double',[xy_dim,time_dim]);
vidconMFf   = netcdf.defVar(ncidMF,'conF','double',[xy_dim,time_dim]);
vidconMFp   = netcdf.defVar(ncidMF,'conP','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidconMPzm  = netcdf.defVar(ncidMP,'conZm','double',[xy_dim,time_dim]);
vidconMPzl  = netcdf.defVar(ncidMP,'conZl','double',[xy_dim,time_dim]);
vidconMPf   = netcdf.defVar(ncidMP,'conF','double',[xy_dim,time_dim]);
vidconMPp   = netcdf.defVar(ncidMP,'conP','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidconMD    = netcdf.defVar(ncidMD,'conB','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidconLPf    = netcdf.defVar(ncidLP,'conF','double',[xy_dim,time_dim]);
vidconLPp    = netcdf.defVar(ncidLP,'conP','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidconLDf    = netcdf.defVar(ncidLD,'conF','double',[xy_dim,time_dim]);
vidconLDp    = netcdf.defVar(ncidLD,'conP','double',[xy_dim,time_dim]);
vidconLDd    = netcdf.defVar(ncidLD,'conD','double',[xy_dim,time_dim]);
vidconLDb    = netcdf.defVar(ncidLD,'conB','double',[xy_dim,time_dim]);
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
        
        S_Sml_f_con(:,DY)   = Sml_f.I;
        S_Sml_p_con(:,DY)   = Sml_p.I;
        S_Sml_d_con(:,DY)   = Sml_d.I;
        S_Med_f_conZm(:,DY) = Med_f.con_zm;
        S_Med_f_conZl(:,DY) = Med_f.con_zl;
        S_Med_f_conF(:,DY)  = Med_f.con_f;
        S_Med_f_conP(:,DY)  = Med_f.con_p;
        S_Med_p_conZm(:,DY) = Med_p.con_zm;
        S_Med_p_conZl(:,DY) = Med_p.con_zl;
        S_Med_p_conF(:,DY)  = Med_p.con_f;
        S_Med_p_conP(:,DY)  = Med_p.con_p;
        S_Med_d_con(:,DY)   = Med_d.I;
        S_Lrg_p_conF(:,DY)  = Lrg_p.con_f;
        S_Lrg_p_conP(:,DY)  = Lrg_p.con_p;
        S_Lrg_d_conF(:,DY)  = Lrg_d.con_f;
        S_Lrg_d_conP(:,DY)  = Lrg_d.con_p;
        S_Lrg_d_conD(:,DY)  = Lrg_d.con_d;
        S_Lrg_d_conB(:,DY)  = Lrg_d.con_be;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidconSF,[0 MNT-1],[NX 1],mean(S_Sml_f_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidconSP,[0 MNT-1],[NX 1],mean(S_Sml_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidconSD,[0 MNT-1],[NX 1],mean(S_Sml_d_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidconMFzm,[0 MNT-1],[NX 1],mean(S_Med_f_conF(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidconMFzl,[0 MNT-1],[NX 1],mean(S_Med_f_conP(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidconMFf,[0 MNT-1],[NX 1],mean(S_Med_f_conF(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidconMFp,[0 MNT-1],[NX 1],mean(S_Med_f_conP(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidconMPzm,[0 MNT-1],[NX 1],mean(S_Med_p_conF(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidconMPzl,[0 MNT-1],[NX 1],mean(S_Med_p_conP(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidconMPf,[0 MNT-1],[NX 1],mean(S_Med_p_conF(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidconMPp,[0 MNT-1],[NX 1],mean(S_Med_p_conP(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidconMD,[0 MNT-1],[NX 1],mean(S_Med_d_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidconLPf,[0 MNT-1],[NX 1],mean(S_Lrg_p_conF(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidconLPp,[0 MNT-1],[NX 1],mean(S_Lrg_p_conP(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidconLDf,[0 MNT-1],[NX 1],mean(S_Lrg_d_conF(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidconLDp,[0 MNT-1],[NX 1],mean(S_Lrg_d_conP(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidconLDd,[0 MNT-1],[NX 1],mean(S_Lrg_d_conD(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidconLDb,[0 MNT-1],[NX 1],mean(S_Lrg_d_conB(:,a(i):b(i)),2));
        
                
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
