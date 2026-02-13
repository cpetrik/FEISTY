% FEISTY output at all locations
% Hist obsclim pristine 1 degree
% Early periods (1841-1960) can have 1841 as the reference year 
% Historic experimental period uses 1901 as the reference year
% Get 1990 for vertical IC file
% Divide fish biomass into 35 vert levels

clear 
close all

load('/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/FishMIP/mfiles_FishMIP/Phase3a/FishMIP_phase3a_exper_times.mat')
time_long_name = 'time';
time_standard_name = 'time';
time_units = 'days since 1901-1-1 00:00:00'; %hist_units;
time_axis = 'T';
calendar = '365_day';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/FishMIP/GFDL_mom6_cobalt2/' cfile '/OneDeg/'];
spath=['/Volumes/petrik-lab/Feisty/NC/MOM6-1D/Global/offline_feisty/'];

%% Benthic material
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

allBB = biomass;
clear biomass

%% find 1990 start index
day90 = 89*365;
id = find(hist_time >= day90);
yid = id(1);

[nid,nt] = size(allBB);

%% SP
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
[nid,nt] = size(biomass);

SP.bio = biomass;
clear biomass

% MP
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass

% LP
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass


% SD
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass 

% MD
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass

% LD
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass


% SF 
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass 

% MF
ncid = netcdf.open([fpath 'Hist_obsclim_pristine_empHP_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,yid-1],[nid 1]);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass


%% Divide up into vertical levels
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%1-D
load([vpath 'grid_OM4_05_COBALTv3.mat'],'z_l','z_l_long_name','z_l_units');
load([vpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'],'thkcello')
thkcello = squeeze(thkcello(:,:,:,1));
thkcello = reshape(thkcello,720*576,35);
thk = mean(thkcello,1,'omitnan');

%% FIsh
SFz = zeros(nid,35);
MFz = zeros(nid,35);
Bz  = zeros(nid,35);
SPz = zeros(nid,35);
MPz = zeros(nid,35);
LPz = zeros(nid,35);
SDz = zeros(nid,35);
MDz = zeros(nid,35);
LDz = zeros(nid,35);

%% Put fish in top 200m
zid = find(z_l<=200);

thkmat = repmat(thk(zid),nid,1);

SFbio = repmat(SF.bio,1,10);
MFbio = repmat(MF.bio,1,10);
SPbio = repmat(SP.bio,1,10);
MPbio = repmat(MP.bio,1,10);
LPbio = repmat(LP.bio,1,10);
SDbio = repmat(SD.bio,1,10);

SFz(:,zid) = SFbio./(10*thkmat);
MFz(:,zid) = MFbio./(10*thkmat);
SPz(:,zid) = SPbio./(10*thkmat);
MPz(:,zid) = MPbio./(10*thkmat);
LPz(:,zid) = LPbio./(10*thkmat);
SDz(:,zid) = SDbio./(10*thkmat);

MDz(:,35) = MD.bio;
LDz(:,35) = LD.bio;
Bz(:,35) = allBB(:,yid);


%% ESM
epath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';
load([epath 'Data_grid_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'],'GRD');
load([epath 'gridspec_gfdl-mom6-cobalt2_obsclim_deptho_onedeg.mat'],'LAT','LON');

[ni,nj] = size(LAT);

lat = LAT(1,:);
lon = LON(:,1);

%% Reshape to lat,lon,yr
Sf_B = 1.0e20*ones(ni,nj,35);
Mf_B = Sf_B;
Sp_B = Sf_B;
Mp_B = Sf_B;
Lp_B = Sf_B;
Sd_B = Sf_B;
Md_B_btf = Sf_B;
Ld_B_btf = Sf_B;
BE_B_btf = Sf_B;

for y=1:10
    gsf = 1.0e20*ones(ni,nj);
    gsf(GRD.ID) = SFz(:,y);
    Sf_B(:,:,y) = gsf;
    
    gmf = 1.0e20*ones(ni,nj);
    gmf(GRD.ID) = MFz(:,y);
    Mf_B(:,:,y) = gmf;
    
    gsp = 1.0e20*ones(ni,nj);
    gsp(GRD.ID) = SPz(:,y);
    Sp_B(:,:,y) = gsp;
    
    gmp = 1.0e20*ones(ni,nj);
    gmp(GRD.ID) = MPz(:,y);
    Mp_B(:,:,y) = gmp;

    glp = 1.0e20*ones(ni,nj);
    glp(GRD.ID) = LPz(:,y);
    Lp_B(:,:,y) = glp;

    gsd = 1.0e20*ones(ni,nj);
    gsd(GRD.ID) = SDz(:,y);
    Sd_B(:,:,y) = gsd;
end


gmd = 1.0e20*ones(ni,nj);
gmd(GRD.ID) = MDz(:,35);
Md_B_btf(:,:,35) = gmd;

gld = 1.0e20*ones(ni,nj);
gld(GRD.ID) = LDz(:,35);
Ld_B_btf(:,:,35) = gld;

gb = 1.0e20*ones(ni,nj);
gb(GRD.ID) = Bz(:,35);
BE_B_btf(:,:,35) = gb;



%% Quick look
pb = Mf_B(:,:,1);
db = Lp_B(:,:,3);
fb = Ld_B_btf(:,:,35);

figure(1)
pcolor(log10(pb))
shading flat
colormap('jet')
colorbar
clim([-3 -1])
title('MF')

figure(2)
pcolor(log10(db))
shading flat
colormap('jet')
colorbar
clim([-3 -1])
title('LP')

figure(3)
pcolor(log10(fb))
shading flat
colormap('jet')
colorbar
clim([-2 2])
title('LD')

%% ========================== netcdf info ===============================

close all

%% Setup netcdf path to store to
fname1 = 'Fish_ICs_from_offline_COBALT1990.nc';

file_tpb = [spath fname1];

[ni,nj,nt] = size(Sf_B);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

%% tpb
ncidSB = netcdf.create(file_tpb,"CLOBBER");

lon_dim = netcdf.defDim(ncidSB,'lon',ni);
lat_dim = netcdf.defDim(ncidSB,'lat',nj);
time_dim = netcdf.defDim(ncidSB,'depth',35);

vidlat = netcdf.defVar(ncidSB,'lat','double',lat_dim);
netcdf.putAtt(ncidSB,vidlat,'long_name','latitude');
netcdf.putAtt(ncidSB,vidlat,'standard_name','lat');
netcdf.putAtt(ncidSB,vidlat,'units','degrees_north');
netcdf.putAtt(ncidSB,vidlat,'axis','Y');

vidlon = netcdf.defVar(ncidSB,'lon','double',lon_dim);
netcdf.putAtt(ncidSB,vidlon,'long_name','longitude');
netcdf.putAtt(ncidSB,vidlon,'standard_name','lon');
netcdf.putAtt(ncidSB,vidlon,'units','degrees_east' );
netcdf.putAtt(ncidSB,vidlon,'axis','X');

vidtSB = netcdf.defVar(ncidSB,'z_l','double',time_dim);
netcdf.putAtt(ncidSB,vidtSB,'long_name',z_l_long_name);
netcdf.putAtt(ncidSB,vidtSB,'units',z_l_units);
netcdf.putAtt(ncidSB,vidtSB,'axis','z_l');

%small
vidbioSF = netcdf.defVar(ncidSB,'SF_B','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSF,'long_name','Small Forage Biomass Density');
netcdf.putAtt(ncidSB,vidbioSF,'units','gWW m-3' );
netcdf.defVarFill(ncidSB,vidbioSF,false,1.0e20);

vidbioSP = netcdf.defVar(ncidSB,'SP_B','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSP,'long_name','Small Lg Pelagic Biomass Density');
netcdf.putAtt(ncidSB,vidbioSP,'units','gWW m-3' );
netcdf.defVarFill(ncidSB,vidbioSP,false,1.0e20);

vidbioSD = netcdf.defVar(ncidSB,'SD_B','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioSD,'long_name','Small Demersal Biomass Density');
netcdf.putAtt(ncidSB,vidbioSD,'units','gWW m-3' );
netcdf.defVarFill(ncidSB,vidbioSD,false,1.0e20);

%medium
vidbioMF = netcdf.defVar(ncidSB,'MF_B','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioMF,'long_name','Medium Forage Biomass Density');
netcdf.putAtt(ncidSB,vidbioMF,'units','gWW m-3' );
netcdf.defVarFill(ncidSB,vidbioMF,false,1.0e20);

vidbioMP = netcdf.defVar(ncidSB,'MP_B','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioMP,'long_name','Medium Lg Pelagic Biomass Density');
netcdf.putAtt(ncidSB,vidbioMP,'units','gWW m-3' );
netcdf.defVarFill(ncidSB,vidbioMP,false,1.0e20);

vidbioMD = netcdf.defVar(ncidSB,'MD_B_btf','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioMD,'long_name','Medium Demersal Biomass Density');
netcdf.putAtt(ncidSB,vidbioMD,'units','gWW m-2' );
netcdf.defVarFill(ncidSB,vidbioMD,false,1.0e20);

%large
vidbioLP = netcdf.defVar(ncidSB,'LP_B','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioLP,'long_name','Large Lg Pelagic Biomass Density');
netcdf.putAtt(ncidSB,vidbioLP,'units','gWW m-3' );
netcdf.defVarFill(ncidSB,vidbioLP,false,1.0e20);

vidbioLD = netcdf.defVar(ncidSB,'LD_B_btf','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioLD,'long_name','Large Demersal Biomass Density');
netcdf.putAtt(ncidSB,vidbioLD,'units','gWW m-2' );
netcdf.defVarFill(ncidSB,vidbioLD,false,1.0e20);

%Bent
vidbioB = netcdf.defVar(ncidSB,'BE_B_btf','double',[lon_dim,lat_dim,time_dim]);
netcdf.putAtt(ncidSB,vidbioB,'long_name','Total Pelagic Biomass Density');
netcdf.putAtt(ncidSB,vidbioB,'units','gWW m-2' );
netcdf.defVarFill(ncidSB,vidbioB,false,1.0e20);


varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncidSB,varid,'creation_date',datestr(now));
netcdf.putAtt(ncidSB,varid,'_FillValue',1.00e20);
netcdf.putAtt(ncidSB,varid,'contact','C. Petrik');
netcdf.putAtt(ncidSB,varid,'institution','UC San Diego');
netcdf.putAtt(ncidSB,varid,'wet weight:C ratio','9:1');

netcdf.endDef(ncidSB);

netcdf.putVar(ncidSB,vidlat,lat);
netcdf.putVar(ncidSB,vidlon,lon);
netcdf.putVar(ncidSB,vidbioSF,Sf_B);
netcdf.putVar(ncidSB,vidbioMF,Mf_B);
netcdf.putVar(ncidSB,vidbioSP,Sp_B);
netcdf.putVar(ncidSB,vidbioMP,Mp_B);
netcdf.putVar(ncidSB,vidbioLP,Lp_B);
netcdf.putVar(ncidSB,vidbioSD,Sd_B);
netcdf.putVar(ncidSB,vidbioMD,Md_B_btf);
netcdf.putVar(ncidSB,vidbioLD,Ld_B_btf);
netcdf.putVar(ncidSB,vidbioB,BE_B_btf);
netcdf.putVar(ncidSB,vidtSB,z_l);

netcdf.close(ncidSB);

%%
ncdisp(file_tpb)



