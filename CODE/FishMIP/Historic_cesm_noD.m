%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Historic_cesm_noD()

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
frate = 0.0;
dfrate = frate/365.0;

%! Choose parameters from other models of my own combo
%1=Kiorboe&Hirst, 2=Hartvig, 3=mizer, 4=JC15, NA=mine
cfn=nan;
efn=nan;
mfn=nan;

%! Make core parameters/constants (global)
make_parameters_noD() % make core parameters/constants

%! Grid
load('/Volumes/FEISTY/Fish-MIP/CESM/Data_grid_cesm.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;

%! How long to run the model
YEARS = length(1850:2005);
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
[fname,simname] = sub_fname_hist_cesm_noD(frate);

%! Storage variables
S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);

% S_Sml_f_prod = zeros(NX,DAYS);
% S_Sml_p_prod = zeros(NX,DAYS);
% S_Med_f_prod = zeros(NX,DAYS);
% S_Med_p_prod = zeros(NX,DAYS);
% S_Lrg_p_prod = zeros(NX,DAYS);


%% ! Initialize
init_sim = simname;
load(['/Volumes/FEISTY/NC/FishMIP/CESM1-BEC/',init_sim '/Last_mo_spinup_' init_sim '.mat']);
[Sml_f,Sml_p,Med_f,Med_p,Lrg_p] = sub_init_fish_hist_noD(ID,DAYS,Sml_f,Sml_p,Med_f,Med_p,Lrg_p);
ENVR = sub_init_env(ID);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_sml_f.nc'];
file_sml_p = [fname,'_sml_p.nc'];
file_med_f = [fname,'_med_f.nc'];
file_med_p = [fname,'_med_p.nc'];
file_lrg_p = [fname,'_lrg_p.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
% vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
% vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
% vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
% vidfishMF   = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
% vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
% vidfishMP   = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
% vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
% vidfishLP   = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
vidT        = netcdf.defVar(ncidLP,'time','double',time_dim);
netcdf.endDef(ncidLP);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    %! Load a year's CESM data
    ti = num2str(YR+1849);
    load(['/Volumes/FEISTY/Fish-MIP/CESM/Hist/Data_cesm_hist_',ti,'.mat'],'CESM');
    %ZOOP OFF BY 1E3
    CESM.Zm = CESM.Zm * 1e3;
    CESM.Zl = CESM.Zl * 1e3;

    for DAY = 1:DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ',num2str(mod(DY,365))]
        [Sml_f,Sml_p,Med_f,Med_p,Lrg_p,ENVR] = ...
            sub_futbio_cesm_noD(ID,DY,CESM,ENVR,Sml_f,Sml_p,...
            Med_f,Med_p,Lrg_p,dfrate);

        %! Store
        S_Sml_f(:,DY) = Sml_f.bio;
        S_Sml_p(:,DY) = Sml_p.bio;
        S_Med_f(:,DY) = Med_f.bio;
        S_Med_p(:,DY) = Med_p.bio;
        S_Lrg_p(:,DY) = Lrg_p.bio;

%         S_Sml_f_prod(:,DY) = Sml_f.prod;
%         S_Sml_p_prod(:,DY) = Sml_p.prod;
%         S_Med_f_prod(:,DY) = Med_f.prod;
%         S_Med_p_prod(:,DY) = Med_p.prod;
%         S_Lrg_p_prod(:,DY) = Lrg_p.prod;

        % S_Med_f_fish(:,DY) = Med_f.caught;
        % S_Med_p_fish(:,DY) = Med_p.caught;
        % S_Lrg_p_fish(:,DY) = Lrg_p.caught;

    end %Days

    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker

        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidbioSF,[0 MNT-1],[NX 1],mean(S_Sml_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidbioSP,[0 MNT-1],[NX 1],mean(S_Sml_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidbioMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidbioLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidT,MNT-1,1,MNT);

%         netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(S_Sml_f_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(S_Sml_p_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(S_Med_f_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(S_Med_p_prod(:,a(i):b(i)),2));
%         netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_prod(:,a(i):b(i)),2));

        % netcdf.putVar(ncidMF,vidfishMF,[0 MNT-1],[NX 1],mean(S_Med_f_fish(:,a(i):b(i)),2));
        % netcdf.putVar(ncidMP,vidfishMP,[0 MNT-1],[NX 1],mean(S_Med_p_fish(:,a(i):b(i)),2));
        % netcdf.putVar(ncidLP,vidfishLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_fish(:,a(i):b(i)),2));

    end %Monthly mean

end %Years


%! Close save
netcdf.close(ncidSF);
netcdf.close(ncidSP);
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidLP);

end
