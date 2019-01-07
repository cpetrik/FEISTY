%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_con_rec_rep()

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
fname = sub_fname(frate);

%! Storage variables
S_Sml_f_rec = zeros(NX,DAYS);
S_Sml_p_rec = zeros(NX,DAYS);
S_Sml_d_rec = zeros(NX,DAYS);
S_Med_f_rec = zeros(NX,DAYS);
S_Med_p_rec = zeros(NX,DAYS);
S_Med_d_rec = zeros(NX,DAYS);
S_Lrg_p_rec = zeros(NX,DAYS);
S_Lrg_d_rec = zeros(NX,DAYS);

S_Sml_f_con = zeros(NX,DAYS);
S_Sml_p_con = zeros(NX,DAYS);
S_Sml_d_con = zeros(NX,DAYS);
S_Med_f_con = zeros(NX,DAYS);
S_Med_p_con = zeros(NX,DAYS);
S_Med_d_con = zeros(NX,DAYS);
S_Lrg_p_con = zeros(NX,DAYS);
S_Lrg_d_con = zeros(NX,DAYS);

S_Med_f_rep = zeros(NX,DAYS);
S_Lrg_p_rep = zeros(NX,DAYS);
S_Lrg_d_rep = zeros(NX,DAYS);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
Med_d.td(1:NX) = 0.0;
Lrg_d.td(1:NX) = 0.0;
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_con_rec_rep_sml_f.nc'];
file_sml_p = [fname,'_con_rec_rep_sml_p.nc'];
file_sml_d = [fname,'_con_rec_rep_sml_d.nc'];
file_med_f = [fname,'_con_rec_rep_med_f.nc'];
file_med_p = [fname,'_con_rec_rep_med_p.nc'];
file_med_d = [fname,'_con_rec_rep_med_d.nc'];
file_lrg_p = [fname,'_con_rec_rep_lrg_p.nc'];
file_lrg_d = [fname,'_con_rec_rep_lrg_d.nc'];


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
vidrecSF    = netcdf.defVar(ncidSF,'rec','double',[xy_dim,time_dim]);
vidconSF    = netcdf.defVar(ncidSF,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidrecSP    = netcdf.defVar(ncidSP,'rec','double',[xy_dim,time_dim]);
vidconSP    = netcdf.defVar(ncidSP,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidrecSD    = netcdf.defVar(ncidSD,'rec','double',[xy_dim,time_dim]);
vidconSD    = netcdf.defVar(ncidSD,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidrecMF    = netcdf.defVar(ncidMF,'rec','double',[xy_dim,time_dim]);
vidconMF    = netcdf.defVar(ncidMF,'con','double',[xy_dim,time_dim]);
vidrepMF    = netcdf.defVar(ncidMF,'rep','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidrecMP    = netcdf.defVar(ncidMP,'rec','double',[xy_dim,time_dim]);
vidconMP    = netcdf.defVar(ncidMP,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidrecMD    = netcdf.defVar(ncidMD,'rec','double',[xy_dim,time_dim]);
vidconMD    = netcdf.defVar(ncidMD,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidrecLP    = netcdf.defVar(ncidLP,'rec','double',[xy_dim,time_dim]);
vidconLP    = netcdf.defVar(ncidLP,'con','double',[xy_dim,time_dim]);
vidrepLP    = netcdf.defVar(ncidLP,'rep','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidrecLD    = netcdf.defVar(ncidLD,'rec','double',[xy_dim,time_dim]);
vidconLD    = netcdf.defVar(ncidLD,'con','double',[xy_dim,time_dim]);
vidrepLD    = netcdf.defVar(ncidLD,'rep','double',[xy_dim,time_dim]);
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
            sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate);
        
        S_Sml_f_rec(:,DY) = Sml_f.rec;
        S_Sml_p_rec(:,DY) = Sml_p.rec;
        S_Sml_d_rec(:,DY) = Sml_d.rec;
        S_Med_f_rec(:,DY) = Med_f.rec;
        S_Med_p_rec(:,DY) = Med_p.rec;
        S_Med_d_rec(:,DY) = Med_d.rec;
        S_Lrg_p_rec(:,DY) = Lrg_p.rec;
        S_Lrg_d_rec(:,DY) = Lrg_d.rec;
        
        S_Sml_f_con(:,DY) = Sml_f.I;
        S_Sml_p_con(:,DY) = Sml_p.I;
        S_Sml_d_con(:,DY) = Sml_d.I;
        S_Med_f_con(:,DY) = Med_f.I;
        S_Med_p_con(:,DY) = Med_p.I;
        S_Med_d_con(:,DY) = Med_d.I;
        S_Lrg_p_con(:,DY) = Lrg_p.I;
        S_Lrg_d_con(:,DY) = Lrg_d.I;
        
        S_Med_f_rep(:,DY) = Med_f.rep;
        S_Lrg_p_rep(:,DY) = Lrg_p.rep;
        S_Lrg_d_rep(:,DY) = Lrg_d.rep;
     
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidMF,vidrepMF,[0 MNT-1],[NX 1],mean(S_Med_f_rep(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidrepLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rep(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidrepLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rep(:,a(i):b(i)),2));        

        netcdf.putVar(ncidSF,vidrecSF,[0 MNT-1],[NX 1],mean(S_Sml_f_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidrecSP,[0 MNT-1],[NX 1],mean(S_Sml_p_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidrecSD,[0 MNT-1],[NX 1],mean(S_Sml_d_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidrecMF,[0 MNT-1],[NX 1],mean(S_Med_f_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidrecMP,[0 MNT-1],[NX 1],mean(S_Med_p_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidrecMD,[0 MNT-1],[NX 1],mean(S_Med_d_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidrecLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rec(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidrecLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rec(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidSF,vidconSF,[0 MNT-1],[NX 1],mean(S_Sml_f_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidconSP,[0 MNT-1],[NX 1],mean(S_Sml_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidconSD,[0 MNT-1],[NX 1],mean(S_Sml_d_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidconMF,[0 MNT-1],[NX 1],mean(S_Med_f_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidconMP,[0 MNT-1],[NX 1],mean(S_Med_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidconMD,[0 MNT-1],[NX 1],mean(S_Med_d_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidconLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidconLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_con(:,a(i):b(i)),2));
        
                
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
