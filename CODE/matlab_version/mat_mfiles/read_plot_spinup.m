% POEM output at all locations

clear all
close all

cfile = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0100';

fpath=['/Volumes/GFDL/NC/Matlab_big_size/' cfile '/'];

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/';
ppath = [pp cfile '/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

%% SP
ncid = netcdf.open([fpath 'Spinup_pristine_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
SP.bio = biomass;
SP.clev = clev;
SP.con = con;
SP.die = die;
SP.gamma = gamma;
SP.nu = nu;
SP.prod = prod;
SP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SF
ncid = netcdf.open([fpath 'Spinup_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.clev = clev;
SF.con = con;
SF.die = die;
SF.gamma = gamma;
SF.nu = nu;
SF.prod = prod;
SF.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SD
ncid = netcdf.open([fpath 'Spinup_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.clev = clev;
SD.con = con;
SD.die = die;
SD.gamma = gamma;
SD.nu = nu;
SD.prod = prod;
SD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MP
ncid = netcdf.open([fpath 'Spinup_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.clev = clev;
MP.con = con;
MP.die = die;
MP.gamma = gamma;
MP.nu = nu;
MP.prod = prod;
MP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MF
ncid = netcdf.open([fpath 'Spinup_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.clev = clev;
MF.con = con;
% MF.DD = DD;
MF.die = die;
% MF.egg = egg;
MF.gamma = gamma;
MF.nu = nu;
MF.prod = prod;
MF.rec = rec;
MF.rep = rep;
% MF.S = S;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MD
ncid = netcdf.open([fpath 'Spinup_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.clev = clev;
MD.con = con;
MD.die = die;
MD.gamma = gamma;
MD.nu = nu;
MD.prod = prod;
MD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LP
ncid = netcdf.open([fpath 'Spinup_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.clev = clev;
LP.con = con;
% LP.DD = DD;
LP.die = die;
% LP.egg = egg;
LP.gamma = gamma;
LP.nu = nu;
LP.prod = prod;
LP.rec = rec;
LP.rep = rep;
% LP.S = S;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LD
ncid = netcdf.open([fpath 'Spinup_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.clev = clev;
LD.con = con;
% LD.DD = DD;
LD.die = die;
% LD.egg = egg;
LD.gamma = gamma;
LD.nu = nu;
LD.prod = prod;
LD.rec = rec;
LD.rep = rep;
% LD.S = S;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% Benthic material
ncid = netcdf.open([fpath 'Spinup_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

BENT.bio = biomass;
clear biomass 

%% Take means
% Last year
lyr=time((end-12+1):end);

sp_mlev=mean(SP.clev(:,lyr),2);
sf_mlev=mean(SF.clev(:,lyr),2);
sd_mlev=mean(SD.clev(:,lyr),2);
mp_mlev=mean(MP.clev(:,lyr),2);
mf_mlev=mean(MF.clev(:,lyr),2);
md_mlev=mean(MD.clev(:,lyr),2);
lp_mlev=mean(LP.clev(:,lyr),2);
ld_mlev=mean(LD.clev(:,lyr),2);

sp_gge=mean(SP.nu(:,lyr),2) ./ mean(SP.con(:,lyr),2);
sf_gge=mean(SF.nu(:,lyr),2) ./ mean(SF.con(:,lyr),2);
sd_gge=mean(SD.nu(:,lyr),2) ./ mean(SD.con(:,lyr),2);
mp_gge=mean(MP.nu(:,lyr),2) ./ mean(MP.con(:,lyr),2);
mf_gge=mean(MF.nu(:,lyr),2) ./ mean(MF.con(:,lyr),2);
md_gge=mean(MD.nu(:,lyr),2) ./ mean(MD.con(:,lyr),2);
lp_gge=mean(LP.nu(:,lyr),2) ./ mean(LP.con(:,lyr),2);
ld_gge=mean(LD.nu(:,lyr),2) ./ mean(LD.con(:,lyr),2);


%% Plots in space
[ni,nj]=size(geolon_t);

% Feeding level
Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);

Zsf(grid(:,1))=sf_mlev;
Zsp(grid(:,1))=sp_mlev;
Zsd(grid(:,1))=sd_mlev;
Zmf(grid(:,1))=mf_mlev;
Zmp(grid(:,1))=mp_mlev;
Zmd(grid(:,1))=md_mlev;
Zlp(grid(:,1))=lp_mlev;
Zld(grid(:,1))=ld_mlev;

% sp
figure(1)
surf(geolon_t,geolat_t,Zsp); view(2); hold on;
shading flat
title('mean Larval P feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_SP.png'])

% sf
figure(2)
surf(geolon_t,geolat_t,Zsf); view(2); hold on;
shading flat
title('mean Larval F feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_SF.png'])

% sd
figure(3)
surf(geolon_t,geolat_t,Zsd); view(2); hold on;
shading flat
title('mean Larval D feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_SD.png'])

% mp
figure(4)
surf(geolon_t,geolat_t,Zmp); view(2); hold on;
shading flat
title('mean Juvenile P feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_MP.png'])

% mf
figure(5)
surf(geolon_t,geolat_t,Zmf); view(2); hold on;
shading flat
title('mean Adult F feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_MF.png'])

% md
figure(6)
surf(geolon_t,geolat_t,Zmd); view(2); hold on;
shading flat
title('mean Juvenile D feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_MD.png'])

% lp
figure(7)
surf(geolon_t,geolat_t,Zlp); view(2); hold on;
shading flat
title('mean Adult P feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_LP.png'])

% ld
figure(8)
surf(geolon_t,geolat_t,Zld); view(2); hold on;
shading flat
title('mean Adult D feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_flev_LD.png'])

%% Gross growth efficiency
Gsf=NaN*ones(ni,nj);
Gsp=NaN*ones(ni,nj);
Gsd=NaN*ones(ni,nj);
Gmf=NaN*ones(ni,nj);
Gmp=NaN*ones(ni,nj);
Gmd=NaN*ones(ni,nj);
Glp=NaN*ones(ni,nj);
Gld=NaN*ones(ni,nj);

Gsf(grid(:,1))=sf_gge;
Gsp(grid(:,1))=sp_gge;
Gsd(grid(:,1))=sd_gge;
Gmf(grid(:,1))=mf_gge;
Gmp(grid(:,1))=mp_gge;
Gmd(grid(:,1))=md_gge;
Glp(grid(:,1))=lp_gge;
Gld(grid(:,1))=ld_gge;

% sp
figure(11)
surf(geolon_t,geolat_t,Gsp); view(2); hold on;
shading flat
title('mean Larval P gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_SP.png'])

% sf
figure(12)
surf(geolon_t,geolat_t,Gsf); view(2); hold on;
shading flat
title('mean Larval F gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_SF.png'])

% sd
figure(13)
surf(geolon_t,geolat_t,Gsd); view(2); hold on;
shading flat
title('mean Larval D gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_SD.png'])

% mp
figure(14)
surf(geolon_t,geolat_t,Gmp); view(2); hold on;
shading flat
title('mean Juvenile P gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_MP.png'])

% mf
figure(15)
surf(geolon_t,geolat_t,Gmf); view(2); hold on;
shading flat
title('mean Adult F gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_MF.png'])

% md
figure(16)
surf(geolon_t,geolat_t,Gmd); view(2); hold on;
shading flat
title('mean Juvenile D gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_MD.png'])

%% lp
figure(17)
surf(geolon_t,geolat_t,Glp); view(2); hold on;
shading flat
title('mean Adult P gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_LP.png'])

% ld
figure(18)
surf(geolon_t,geolat_t,Gld); view(2); hold on;
shading flat
title('mean Adult D gross growth efficiency')
colormap('jet')
colorbar('h')
caxis([-0.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_gge_LD.png'])

