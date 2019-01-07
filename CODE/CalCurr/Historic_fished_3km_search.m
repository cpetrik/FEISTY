%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Historic_fished_3km_search()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet 
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global tstep ni nj

%% %%%%%%%%%%%%% Initialize Model Variables
ad = 0.1:0.1:1.0;
% for m=1:length(ad)
% J = ad(m);
    for n=1:length(ad)
    Sm = ad(n);

%! Set fishing rate
frate = 0.1;
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_parameters() % make core parameters/constants

%! Grid
load('/Volumes/GFDL/NEMURO/3km/Data_grid_3km_hist.mat');
NX = length(GRD.Z);
ID = 1:NX;

%! How long to run the model
YEARS = 23; %1988-2010
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
[fname,simname] = sub_fname_hist3km(frate);

%% %%%%%%%%%%%%% Setup save
%! Dims 
nt = 12*YEARS;

xy_dim      = NX;
time_dim    = nt;

%! Monthly Storage variables
SF.bio    = NaN*ones(xy_dim,time_dim);
SP.bio    = NaN*ones(xy_dim,time_dim);
SD.bio    = NaN*ones(xy_dim,time_dim);
MF.bio    = NaN*ones(xy_dim,time_dim);
MF.yield  = NaN*ones(xy_dim,time_dim);
MP.bio    = NaN*ones(xy_dim,time_dim);
MP.yield  = NaN*ones(xy_dim,time_dim);
MD.bio    = NaN*ones(xy_dim,time_dim);
MD.yield  = NaN*ones(xy_dim,time_dim);
LP.bio    = NaN*ones(xy_dim,time_dim);
LP.yield  = NaN*ones(xy_dim,time_dim);
LD.bio    = NaN*ones(xy_dim,time_dim);
LD.yield  = NaN*ones(xy_dim,time_dim);
Bent.bio  = NaN*ones(xy_dim,time_dim);
tmo       = NaN*ones(time_dim,1);

%! Daily Storage variables
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


%% ! Initialize
%init_sim = simname;
init_sim = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
load(['/Volumes/GFDL/NC/Matlab_new_size/',init_sim '/CalCurr/Last_mo_Spinup3km_All_fish03_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
Med_d.td(1:NX) = 0.0;
Lrg_d.td(1:NX) = 0.0;
ENVR = sub_init_env(ID);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    %! Load a year's NEMURO data
    ti = num2str(YR+1987);
    load(['/Volumes/GFDL/NEMURO/3km/daily3km/Data_3km_',ti,'.mat']);
    
    [num2str(YR)]
    
    for DAY = 1:DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(ID,DY,NEMURO,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
        
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
        Bent.bio(:,MNT) = mean(S_Bent_bio(:,a(i):b(i)),2);
        tmo(MNT) = MNT;
        
        SF.bio(:,MNT) = mean(S_Sml_f(:,a(i):b(i)),2);
        SP.bio(:,MNT) = mean(S_Sml_p(:,a(i):b(i)),2);
        SD.bio(:,MNT) = mean(S_Sml_d(:,a(i):b(i)),2);
        MF.bio(:,MNT) = mean(S_Med_f(:,a(i):b(i)),2);
        MP.bio(:,MNT) = mean(S_Med_p(:,a(i):b(i)),2);
        MD.bio(:,MNT) = mean(S_Med_d(:,a(i):b(i)),2);
        LP.bio(:,MNT) = mean(S_Lrg_p(:,a(i):b(i)),2);
        LD.bio(:,MNT) = mean(S_Lrg_d(:,a(i):b(i)),2);
        
        MF.yield(:,MNT) = mean(S_Med_f_fish(:,a(i):b(i)),2);
        MP.yield(:,MNT) = mean(S_Med_p_fish(:,a(i):b(i)),2);
        MD.yield(:,MNT) = mean(S_Med_d_fish(:,a(i):b(i)),2);
        LP.yield(:,MNT) = mean(S_Lrg_p_fish(:,a(i):b(i)),2);
        LD.yield(:,MNT) = mean(S_Lrg_d_fish(:,a(i):b(i)),2);
        
    end %Monthly mean
    
end %Years


%% %%%%%%%%%%%%%%%%%%%% Time and space means
%! Time
sp_tmean=mean(SP.bio,1);
sf_tmean=mean(SF.bio,1);
sd_tmean=mean(SD.bio,1);
mp_tmean=mean(MP.bio,1);
mf_tmean=mean(MF.bio,1);
md_tmean=mean(MD.bio,1);
lp_tmean=mean(LP.bio,1);
ld_tmean=mean(LD.bio,1);
b_tmean=mean(Bent.bio,1);

mf_tmy=mean(MF.yield,1);
mp_tmy=mean(MP.yield,1);
md_tmy=mean(MD.yield,1);
lp_tmy=mean(LP.yield,1);
ld_tmy=mean(LD.yield,1);

%! All years
sp_mean=mean(SP.bio,2);
sf_mean=mean(SF.bio,2);
sd_mean=mean(SD.bio,2);
mp_mean=mean(MP.bio,2);
mf_mean=mean(MF.bio,2);
md_mean=mean(MD.bio,2);
lp_mean=mean(LP.bio,2);
ld_mean=mean(LD.bio,2);
b_mean=mean(Bent.bio,2);

mf_my=mean(MF.yield,2);
mp_my=mean(MP.yield,2);
md_my=mean(MD.yield,2);
lp_my=mean(LP.yield,2);
ld_my=mean(LD.yield,2);

%! Each year
a = 1:12:nt; % start of each yr
b = 12:12:nt; % end of each yr
mB = NaN*ones(length(b_mean),(nt/12));
mSF = mB;
mSP = mB;
mSD = mB;
mMF = mB;
mMP = mB;
mMD = mB;
mLP = mB;
mLD = mB;
for i = 1:(nt/12)
    %yr = (i+1987);
    %! Put vars of netcdf file
    mB(:,i) = mean(Bent.bio(:,a(i):b(i)),2);
    mSF(:,i) = mean(SF.bio(:,a(i):b(i)),2);
    mSP(:,i) = mean(SP.bio(:,a(i):b(i)),2);
    mSD(:,i) = mean(SD.bio(:,a(i):b(i)),2);
    mMF(:,i) = mean(MF.bio(:,a(i):b(i)),2);
    mMP(:,i) = mean(MP.bio(:,a(i):b(i)),2);
    mMD(:,i) = mean(MD.bio(:,a(i):b(i)),2);
    mLP(:,i) = mean(LP.bio(:,a(i):b(i)),2);
    mLD(:,i) = mean(LD.bio(:,a(i):b(i)),2);
end

%! Each season
a = 1:3:nt; % start of each 3-mo period
b = 3:3:nt; % end of each 3-mo period
sB = NaN*ones(length(b_mean),(nt/12));
sSF = sB;
sSP = sB;
sSD = sB;
sMF = sB;
sMP = sB;
sMD = sB;
sLP = sB;
sLD = sB;
for i = 1:(nt/3)
    sB(:,i) = mean(Bent.bio(:,a(i):b(i)),2);
    sSF(:,i) = mean(SF.bio(:,a(i):b(i)),2);
    sSP(:,i) = mean(SP.bio(:,a(i):b(i)),2);
    sSD(:,i) = mean(SD.bio(:,a(i):b(i)),2);
    sMF(:,i) = mean(MF.bio(:,a(i):b(i)),2);
    sMP(:,i) = mean(MP.bio(:,a(i):b(i)),2);
    sMD(:,i) = mean(MD.bio(:,a(i):b(i)),2);
    sLP(:,i) = mean(LP.bio(:,a(i):b(i)),2);
    sLD(:,i) = mean(LD.bio(:,a(i):b(i)),2);
end
    
%! Save
save([fname '_' simname '_Means.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
    'md_tmean','b_tmean','lp_tmean','ld_tmean','tmo',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
    'mf_my','mp_my','md_my','lp_my','ld_my',...
    'mB','mSF','mSP','mSD','mMF','mMP','mMD','mLP','mLD',...
    'sB','sSF','sSP','sSD','sMF','sMP','sMD','sLP','sLD');


    end
% end

end
