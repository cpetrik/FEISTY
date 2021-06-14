%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_death_vars()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep K CGRD ni nj MZpref

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
make_parameters_SW()

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

%! Create a directory for output
%fname = sub_fname(frate);
test = 'SWmarth_MZ25';
fname = sub_fname_test(frate,test);

%! Storage variables
S_Sml_f_mort = zeros(NX,DAYS);
S_Sml_p_mort = zeros(NX,DAYS);
S_Sml_d_mort = zeros(NX,DAYS);
S_Med_f_mort = zeros(NX,DAYS);
S_Med_p_mort = zeros(NX,DAYS);
S_Med_d_mort = zeros(NX,DAYS);
S_Lrg_p_mort = zeros(NX,DAYS);
S_Lrg_d_mort = zeros(NX,DAYS);

S_Sml_f_die = zeros(NX,DAYS);
S_Sml_p_die = zeros(NX,DAYS);
S_Sml_d_die = zeros(NX,DAYS);
S_Med_f_die = zeros(NX,DAYS);
S_Med_p_die = zeros(NX,DAYS);
S_Med_d_die = zeros(NX,DAYS);
S_Lrg_p_die = zeros(NX,DAYS);
S_Lrg_d_die = zeros(NX,DAYS);

S_Med_f_fish = zeros(NX,DAYS);
S_Med_p_fish = zeros(NX,DAYS);
S_Med_d_fish = zeros(NX,DAYS);
S_Lrg_p_fish = zeros(NX,DAYS);
S_Lrg_d_fish = zeros(NX,DAYS);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
Med_d.td(1:NX) = 0.0;
Lrg_d.td(1:NX) = 0.0;
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_die_nmort_yield_sml_f.nc'];
file_sml_p = [fname,'_die_nmort_yield_sml_p.nc'];
file_sml_d = [fname,'_die_nmort_yield_sml_d.nc'];
file_med_f = [fname,'_die_nmort_yield_med_f.nc'];
file_med_p = [fname,'_die_nmort_yield_med_p.nc'];
file_med_d = [fname,'_die_nmort_yield_med_d.nc'];
file_lrg_p = [fname,'_die_nmort_yield_lrg_p.nc'];
file_lrg_d = [fname,'_die_nmort_yield_lrg_d.nc'];

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
vidmortSF     = netcdf.defVar(ncidSF,'mort','double',[xy_dim,time_dim]);
viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidmortSP     = netcdf.defVar(ncidSP,'mort','double',[xy_dim,time_dim]);
viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidmortSD     = netcdf.defVar(ncidSD,'mort','double',[xy_dim,time_dim]);
viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidmortMF     = netcdf.defVar(ncidMF,'mort','double',[xy_dim,time_dim]);
viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
vidfishMF   = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidmortMP     = netcdf.defVar(ncidMP,'mort','double',[xy_dim,time_dim]);
viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
vidfishMP   = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidmortMD     = netcdf.defVar(ncidMD,'mort','double',[xy_dim,time_dim]);
viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
vidfishMD   = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidmortLP     = netcdf.defVar(ncidLP,'mort','double',[xy_dim,time_dim]);
viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
vidfishLP   = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidmortLD     = netcdf.defVar(ncidLD,'mort','double',[xy_dim,time_dim]);
viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
vidfishLD   = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);


%% %%%%%%%%%%%%%%%%%%%% Run the Model %%%%%%%%%%%%%%%%%%%%
%! Run model 
MNT=0;
for YR = 1:YEARS % years
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_mzpref(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
        
        S_Sml_f_die(:,DY) = Sml_f.die;
        S_Sml_p_die(:,DY) = Sml_p.die;
        S_Sml_d_die(:,DY) = Sml_d.die;
        S_Med_f_die(:,DY) = Med_f.die;
        S_Med_p_die(:,DY) = Med_p.die;
        S_Med_d_die(:,DY) = Med_d.die;
        S_Lrg_p_die(:,DY) = Lrg_p.die;
        S_Lrg_d_die(:,DY) = Lrg_d.die;

        S_Sml_f_mort(:,DY) = Sml_f.mort;
        S_Sml_p_mort(:,DY) = Sml_p.mort;
        S_Sml_d_mort(:,DY) = Sml_d.mort;
        S_Med_f_mort(:,DY) = Med_f.mort;
        S_Med_p_mort(:,DY) = Med_p.mort;
        S_Med_d_mort(:,DY) = Med_d.mort;
        S_Lrg_p_mort(:,DY) = Lrg_p.mort;
        S_Lrg_d_mort(:,DY) = Lrg_d.mort;
       
        S_Med_f_fish(:,DY) = Med_f.caught;
        S_Med_p_fish(:,DY) = Med_p.caught;
        S_Med_d_fish(:,DY) = Med_d.caught;
        S_Lrg_p_fish(:,DY) = Lrg_p.caught;
        S_Lrg_d_fish(:,DY) = Lrg_d.caught;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidmortSF,[0 MNT-1],[NX 1],mean(S_Sml_f_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidmortSP,[0 MNT-1],[NX 1],mean(S_Sml_p_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidmortSD,[0 MNT-1],[NX 1],mean(S_Sml_d_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidmortMF,[0 MNT-1],[NX 1],mean(S_Med_f_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidmortMP,[0 MNT-1],[NX 1],mean(S_Med_p_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidmortMD,[0 MNT-1],[NX 1],mean(S_Med_d_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidmortLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidmortLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_mort(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidSF,viddieSF,[0 MNT-1],[NX 1],mean(S_Sml_f_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,viddieSP,[0 MNT-1],[NX 1],mean(S_Sml_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,viddieSD,[0 MNT-1],[NX 1],mean(S_Sml_d_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,viddieMF,[0 MNT-1],[NX 1],mean(S_Med_f_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,viddieMP,[0 MNT-1],[NX 1],mean(S_Med_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,viddieMD,[0 MNT-1],[NX 1],mean(S_Med_d_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,viddieLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,viddieLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_die(:,a(i):b(i)),2));
           
        netcdf.putVar(ncidMF,vidfishMF,[0 MNT-1],[NX 1],mean(S_Med_f_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidfishMP,[0 MNT-1],[NX 1],mean(S_Med_p_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidfishMD,[0 MNT-1],[NX 1],mean(S_Med_d_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidfishLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_fish(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidfishLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_fish(:,a(i):b(i)),2));
        
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
