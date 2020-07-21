%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_fished_Zloss()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj kc

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
frate = 0.3; %Fish(F);
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_parameters()

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Grid; choose where and when to run the model
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

%! Create a directory for output
fname = sub_fname(frate);

%! Storage variables
S_Mzoo_frac = zeros(NX,DAYS);
S_Lzoo_frac = zeros(NX,DAYS);
S_Bent_frac = zeros(NX,DAYS);


%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = ...
        sub_init_fish(ID,DAYS);
Med_d.td(1:NX) = 0.0;
Lrg_d.td(1:NX) = 0.0;
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_mzoo  = [fname,'_mzoo.nc'];
file_lzoo  = [fname,'_lzoo.nc'];
file_bent  = [fname,'_bfrac.nc'];

ncidMZ = netcdf.create(file_mzoo,'NC_WRITE');
ncidLZ = netcdf.create(file_lzoo,'NC_WRITE');
ncidB  = netcdf.create(file_bent,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
%oldFormat = netcdf.setDefaultFormat('NC_FORMAT_64BIT');
netcdf.setDefaultFormat('NC_FORMAT_64BIT')

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidMZ,'nid',NX);
time_dim    = netcdf.defDim(ncidMZ,'ntime',nt);
vidfracMZ   = netcdf.defVar(ncidMZ,'fraction','double',[xy_dim,time_dim]);
vidTMZ      = netcdf.defVar(ncidMZ,'time','double',time_dim);
netcdf.endDef(ncidMZ);

xy_dim      = netcdf.defDim(ncidLZ,'nid',NX);
time_dim    = netcdf.defDim(ncidLZ,'ntime',nt);
vidfracLZ    = netcdf.defVar(ncidLZ,'fraction','double',[xy_dim,time_dim]);
vidTLZ      = netcdf.defVar(ncidLZ,'time','double',time_dim);
netcdf.endDef(ncidLZ);

xy_dim     = netcdf.defDim(ncidB,'nid',NX);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'fraction','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);

%% %%%%%%%%%%%%%%%%%%%% Run the Model %%%%%%%%%%%%%%%%%%%%
%! Run model 
MNT=0;
for YR = 1:YEARS % years
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
        
        S_Mzoo_frac(:,DY) = ENVR.fZm;
        S_Lzoo_frac(:,DY) = ENVR.fZl;
        S_Bent_frac(:,DY) = ENVR.fB;
        
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidB,vidbioB,[0 MNT-1],[NX 1],mean(S_Bent_frac(:,a(i):b(i)),2));
        netcdf.putVar(ncidB,vidTB,MNT-1,1,MNT);
        
        netcdf.putVar(ncidMZ,vidfracMZ,[0 MNT-1],[NX 1],mean(S_Mzoo_frac(:,a(i):b(i)),2));
        netcdf.putVar(ncidLZ,vidfracLZ,[0 MNT-1],[NX 1],mean(S_Lzoo_frac(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidMZ,vidTMZ,MNT-1,1,MNT);
        netcdf.putVar(ncidLZ,vidTLZ,MNT-1,1,MNT);
        
    end %Monthly mean
    
end %Years

%! Close save
netcdf.close(ncidMZ);
netcdf.close(ncidLZ);
netcdf.close(ncidB);

end
