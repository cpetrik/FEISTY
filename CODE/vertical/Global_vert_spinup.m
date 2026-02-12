%%%%!! RUN SPINUP FOR WHOLE GLOBE, ALL DEPTHS
function Global_vert_spinup()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;
param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_parameters(param);

%! Grids
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
%vpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
%vpath = '/scratch/cpetrik/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%1-D
load([vpath 'grid_OM4_05_COBALTv3.mat'],'wet','geolon','z_l');

param.depth = 4.0e3; %get from grid eventually
param.hploss = 1; %1=yes cap; 0=no cap

%ocean grid cells
WID = find(wet(:)==1);
NWID = length(WID);
param.WID = WID;
param.NWID = 1:NWID;
param.NW = NWID;
%2Dgrid
[ni,nj] = size(geolon);
%vertical layers
NZID = length(z_l);
param.ZID = 1:NZID;
param.NZ = NZID;
param.NX = NZID;

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
aa = (cumsum(MNTH)+1);
a = [1,aa(1:end-1)]; %start of the month
b = cumsum(MNTH);    %end of the month
Tdays=1:DAYS;

%! Create a directory for output
%eventually change so exper is subfolder within offline_feisty
exper = 'Global_spinup_COBALTv3_halfdeg';
opath = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/Global/offline_feisty/';
%opath = '/project/Feisty/NC/MOM6-1D/Global/offline_feisty/';
[fname,simname,sname] = sub_fname_spin(param,opath,exper);


%! Dims of netcdf file
% Too big for monthly or annual means, do every 5 yrs
nt = YEARS/5;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% %%%%%%%%%%%%% Setup NetCDF save
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

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');
ncidB  = netcdf.create(file_bent,'NC_WRITE');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~15 minutes ... ']
%3-D
xy_dim      = netcdf.defDim(ncidSF,'wid',NWID);
z_dim       = netcdf.defDim(ncidSF,'zid',NZID);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,z_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'wid',NWID);
z_dim       = netcdf.defDim(ncidSP,'zid',NZID);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,z_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'wid',NWID);
z_dim       = netcdf.defDim(ncidSD,'zid',NZID);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,z_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'wid',NWID);
z_dim       = netcdf.defDim(ncidMF,'zid',NZID);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,z_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'wid',NWID);
z_dim       = netcdf.defDim(ncidMP,'zid',NZID);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,z_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidLP,'wid',NWID);
z_dim       = netcdf.defDim(ncidLP,'zid',NZID);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,z_dim,time_dim]);
netcdf.endDef(ncidLP);

%2-D
xy_dim      = netcdf.defDim(ncidMD,'wid',NWID);
%z_dim       = netcdf.defDim(ncidMD,'zid',NZID);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
%vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,z_dim,time_dim]);
vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLD,'wid',NWID);
%z_dim       = netcdf.defDim(ncidLD,'zid',NZID);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
%vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,z_dim,time_dim]);
vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

xy_dim      = netcdf.defDim(ncidB,'wid',NWID);
time_dim   = netcdf.defDim(ncidB,'ntime',nt);
vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
netcdf.endDef(ncidB);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Initial year of forcing

%READ IN YEAR OF INTEREST ONLY FROM NETCDF (only for spinup, read all 5 for hist)
[TEMP_z,TEMP_btm,MZ,LZ,MZloss,LZloss,det_btm] = ncread_global_feisty_forcing_first_year(vpath,ni,nj,NZID);

% THICKNESS
load([vpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'],'thkcello')
thkcello = thkcello(:,:,:,1:12);

%! Loop over ocean grid cells - PARALLELIZE THIS STEP
%parfor W = 1:NWID
for W = 1:4%:NWID

    [m,n] = ind2sub([ni,nj],WID(W)); % spatial index of water cell
    % location of interest
    Tp  = double(squeeze(TEMP_z(m,n,:,:)));
    Zm  = double(squeeze(MZ(m,n,:,:)));
    Zl  = double(squeeze(LZ(m,n,:,:)));
    dZm = double(squeeze(MZloss(m,n,:,:)));
    dZl = double(squeeze(LZloss(m,n,:,:)));
    thk = double(squeeze(thkcello(m,n,:,:)));
    Tb  = double(squeeze(TEMP_btm(m,n,:)));
    det = double(squeeze(det_btm(m,n,:)));

    Tb = Tb';
    det = det';

    %Interpolate monthly forcing to daily
    Time=15:30:DAYS;
    ESM = daily_interp_monthly_means(NZID,Time,Tdays,...
        Tp,Tb,det,Zm,Zl,dZm,dZl,thk);

    %! Initialize
    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(NZID,param);

    %! Iterate forward in time
    yrtick = 0;
    for YR = 1:YEARS % years
        %num2str(YR)

        %! Storage
        D_Sml_f = zeros(NZID,DAYS);
        D_Sml_p = zeros(NZID,DAYS);
        D_Sml_d = zeros(NZID,DAYS);
        D_Med_f = zeros(NZID,DAYS);
        D_Med_p = zeros(NZID,DAYS);
        D_Lrg_p = zeros(NZID,DAYS);

        D_Med_d = zeros(1,DAYS);
        D_Lrg_d = zeros(1,DAYS);
        D_Bent_bio = zeros(1,DAYS);

        % SIMULATE THE DAYS IN EACH MONTH INSTEAD OF WHOLE YEAR?
        for DY = 1:DAYS % days

            %%%! Future time step
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                sub_futbio_glob(DY,ESM,Sml_f,Sml_p,Sml_d,...
                Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

            %! Daily for one year
            D_Sml_f(:,DY) = Sml_f.bio;
            D_Sml_p(:,DY) = Sml_p.bio;
            D_Sml_d(:,DY) = Sml_d.bio;
            D_Med_f(:,DY) = Med_f.bio;
            D_Med_p(:,DY) = Med_p.bio;
            D_Lrg_p(:,DY) = Lrg_p.bio;

            D_Bent_bio(:,DY) = BENT.mass(NZID);
            D_Med_d(:,DY) = Med_d.bio(NZID);
            D_Lrg_d(:,DY) = Lrg_d.bio(NZID);


        end %Days

        %! Calculate annual means and save
        %! Put vars of netcdf file
        if (rem(YR,5)==0) % every 5 yrs
            yrtick = yrtick +1;
            %3-D
            netcdf.putVar(ncidSF,vidbioSF,[W-1,0,yrtick-1],[1,NZID,1],mean(D_Sml_f,2));
            netcdf.putVar(ncidSP,vidbioSP,[W-1,0,yrtick-1],[1,NZID,1],mean(D_Sml_p,2));
            netcdf.putVar(ncidSD,vidbioSD,[W-1,0,yrtick-1],[1,NZID,1],mean(D_Sml_d,2));
            netcdf.putVar(ncidMF,vidbioMF,[W-1,0,yrtick-1],[1,NZID,1],mean(D_Med_f,2));
            netcdf.putVar(ncidMP,vidbioMP,[W-1,0,yrtick-1],[1,NZID,1],mean(D_Med_p,2));
            netcdf.putVar(ncidLP,vidbioLP,[W-1,0,yrtick-1],[1,NZID,1],mean(D_Lrg_p,2));
            %2-D
            netcdf.putVar(ncidMD,vidbioMD,[W-1,yrtick-1],[1,1],mean(D_Med_d,2));
            netcdf.putVar(ncidLD,vidbioLD,[W-1,yrtick-1],[1,1],mean(D_Lrg_d,2));
            netcdf.putVar(ncidB,vidbioB,[W-1,yrtick-1],[1,1],mean(D_Bent_bio,2));
            netcdf.putVar(ncidB,vidTB,yrtick-1,1,YR);
        end % every 5 yrs

    end %Years

end %Grid cells

%%% Save NETCDFS
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


end