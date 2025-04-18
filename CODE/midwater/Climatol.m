%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants (global)
param = [];
param = make_parameters(param);

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Grid; choose where and when to run the model
load('/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;
param.NX = NX;
param.ID = ID;

%! Create a directory for output
Fdir = pwd;
addpath([Fdir '/naming/'])
[fname] = sub_fname_clim(param);

%! Storage variables
S_Bent_bio = zeros(NX,DAYS);
S_Mzoo_frac = zeros(NX,DAYS);
S_Lzoo_frac = zeros(NX,DAYS);

S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Sml_d = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

S_Sml_f_prod = zeros(NX,DAYS);
S_Sml_p_prod = zeros(NX,DAYS);
S_Sml_d_prod = zeros(NX,DAYS);
S_Med_f_prod = zeros(NX,DAYS);
S_Med_p_prod = zeros(NX,DAYS);
S_Med_d_prod = zeros(NX,DAYS);
S_Lrg_p_prod = zeros(NX,DAYS);
S_Lrg_d_prod = zeros(NX,DAYS);

S_Med_f_fish = zeros(NX,DAYS);
S_Med_p_fish = zeros(NX,DAYS);
S_Med_d_fish = zeros(NX,DAYS);
S_Lrg_p_fish = zeros(NX,DAYS);
S_Lrg_d_fish = zeros(NX,DAYS);


%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = ...
        sub_init_fish(ID);

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
file_mzoo  = [fname,'_mzoo.nc'];
file_lzoo  = [fname,'_lzoo.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');
ncidB  = netcdf.create(file_bent,'NC_WRITE');
ncidMZ = netcdf.create(file_mzoo,'NC_WRITE');
ncidLZ = netcdf.create(file_lzoo,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT')

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,time_dim]);
vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
vidfishMF   = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
vidfishMP   = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
vidfishMD   = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
vidfishLP   = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
vidfishLD   = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

xy_dim     = netcdf.defDim(ncidB,'nid',NX);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);

xy_dim      = netcdf.defDim(ncidMZ,'nid',NX);
time_dim    = netcdf.defDim(ncidMZ,'ntime',nt);
vidfracMZ   = netcdf.defVar(ncidMZ,'fraction','double',[xy_dim,time_dim]);
vidTMZ      = netcdf.defVar(ncidMZ,'time','double',time_dim);
netcdf.endDef(ncidMZ);

xy_dim      = netcdf.defDim(ncidLZ,'nid',NX);
time_dim    = netcdf.defDim(ncidLZ,'ntime',nt);
vidfracLZ   = netcdf.defVar(ncidLZ,'fraction','double',[xy_dim,time_dim]);
vidTLZ      = netcdf.defVar(ncidLZ,'time','double',time_dim);
netcdf.endDef(ncidLZ);

%% %%%%%%%%%%%%%%%%%%%% Run the Model %%%%%%%%%%%%%%%%%%%%
%! Run model 
MNT=0;
for YR = 1:YEARS % years
    
    num2str(YR)
    
    for DAY = 1:param.DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(ID,DY,COBALT,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
        
        S_Bent_bio(:,DY) = BENT.mass;
        S_Mzoo_frac(:,DY) = ENVR.fZm;
        S_Lzoo_frac(:,DY) = ENVR.fZl;
        
        S_Sml_f(:,DY) = Sml_f.bio;
        S_Sml_p(:,DY) = Sml_p.bio;
        S_Sml_d(:,DY) = Sml_d.bio;
        S_Med_f(:,DY) = Med_f.bio;
        S_Med_p(:,DY) = Med_p.bio;
        S_Med_d(:,DY) = Med_d.bio;
        S_Lrg_p(:,DY) = Lrg_p.bio;
        S_Lrg_d(:,DY) = Lrg_d.bio;
        
        S_Sml_f_prod(:,DY) = Sml_f.prod;
        S_Sml_p_prod(:,DY) = Sml_p.prod;
        S_Sml_d_prod(:,DY) = Sml_d.prod;
        S_Med_f_prod(:,DY) = Med_f.prod;
        S_Med_p_prod(:,DY) = Med_p.prod;
        S_Med_d_prod(:,DY) = Med_d.prod;
        S_Lrg_p_prod(:,DY) = Lrg_p.prod;
        S_Lrg_d_prod(:,DY) = Lrg_d.prod;
        
        S_Med_f_fish(:,DY) = Med_f.caught;
        S_Med_p_fish(:,DY) = Med_p.caught;
        S_Med_d_fish(:,DY) = Med_d.caught;
        S_Lrg_p_fish(:,DY) = Lrg_p.caught;
        S_Lrg_d_fish(:,DY) = Lrg_d.caught;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH);    % end of the month
    for i = 1:12
        MNT = MNT+1;     % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidB,vidbioB,[0 MNT-1],[NX 1],mean(S_Bent_bio(:,a(i):b(i)),2));
        netcdf.putVar(ncidB,vidTB,MNT-1,1,MNT);
        
        netcdf.putVar(ncidMZ,vidfracMZ,[0 MNT-1],[NX 1],mean(S_Mzoo_frac(:,a(i):b(i)),2));
        netcdf.putVar(ncidMZ,vidTMZ,MNT-1,1,MNT);
        
        netcdf.putVar(ncidLZ,vidfracLZ,[0 MNT-1],[NX 1],mean(S_Lzoo_frac(:,a(i):b(i)),2));
        netcdf.putVar(ncidLZ,vidTLZ,MNT-1,1,MNT);
        
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
        
        netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(S_Sml_f_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(S_Sml_p_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidprodSD,[0 MNT-1],[NX 1],mean(S_Sml_d_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(S_Med_f_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(S_Med_p_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidprodMD,[0 MNT-1],[NX 1],mean(S_Med_d_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_prod(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidprodLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_prod(:,a(i):b(i)),2));
        
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
netcdf.close(ncidMZ);
netcdf.close(ncidLZ);


end
