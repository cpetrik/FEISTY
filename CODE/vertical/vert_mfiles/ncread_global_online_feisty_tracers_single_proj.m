% 0.5 degree global MOM6-COBALTv3-FEISTY online


clear
close all

%fpath = '/Volumes/petrik-lab/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';
fpath = '/project/Feisty/NC/Global_COBALT_FEISTY/cobalt_feisty/';

%gpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';
gpath = '/project/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
%load([gpath 'grid_OM4_05_COBALTv3.mat'],'wet','z_l_units','z_l_long_name','z_l')

%%
%ncdisp([fpath '19900101.ocean_feisty_tracers_z.nc'])

%%
ncid = netcdf.open([fpath '19900101.ocean_feisty_tracers_z.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:2
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
% 3 & 4 are '01_l' & '01_i' = can't start matlab var with 0
for i = 5:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end
netcdf.close(ncid);

Sf_B(Sf_B>1e19) = nan; 
Sp_B(Sp_B>1e19) = nan;
Sd_B(Sd_B>1e19) = nan;
Mf_B(Mf_B>1e19) = nan;
Mp_B(Mp_B>1e19) = nan;
Md_B(Md_B>1e19) = nan;
Lp_B(Lp_B>1e19) = nan;
Ld_B(Ld_B>1e19) = nan;
BE_B(BE_B>1e19) = nan; % >1e15?

%% thkcello
load([gpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'])

thkcello = thkcello(:,:,:,1:12);

%% Vertical means
vSf = mean(Sf_B,1,'omitnan');
vSf = squeeze(mean(vSf,2,'omitnan'));
vSF = squeeze(mean(vSf,2,'omitnan'));

vSp = mean(Sp_B,1,'omitnan');
vSp = squeeze(mean(vSp,2,'omitnan'));
vSP = squeeze(mean(vSp,2,'omitnan'));

vSd = mean(Sd_B,1,'omitnan');
vSd = squeeze(mean(vSd,2,'omitnan'));
vSD = squeeze(mean(vSd,2,'omitnan'));

vMf = mean(Mf_B,1,'omitnan');
vMf = squeeze(mean(vMf,2,'omitnan'));
vMF = squeeze(mean(vMf,2,'omitnan'));

vMp = mean(Mp_B,1,'omitnan');
vMp = squeeze(mean(vMp,2,'omitnan'));
vMP = squeeze(mean(vMp,2,'omitnan'));

vLp = mean(Lp_B,1,'omitnan');
vLp = squeeze(mean(vLp,2,'omitnan'));
vLP = squeeze(mean(vLp,2,'omitnan'));

% sum btm vars?
vMd = mean(Md_B,1,'omitnan');
vMd = squeeze(mean(vMd,2,'omitnan'));
vMD = squeeze(mean(vMd,2,'omitnan'));
%vMD = squeeze(sum(vMd,2,'omitnan'));

vLd = mean(Ld_B,1,'omitnan');
vLd = squeeze(mean(vLd,2,'omitnan'));
vLD = squeeze(mean(vLd,2,'omitnan'));

vBe = mean(BE_B,1,'omitnan');
vBe = squeeze(mean(vBe,2,'omitnan'));
vBE = squeeze(mean(vBe,2,'omitnan'));

%% vert sums
iSf = squeeze(sum((Sf_B.*thkcello),3,'omitnan'));
iMf = squeeze(sum((Mf_B.*thkcello),3,'omitnan'));
iSp = squeeze(sum((Sp_B.*thkcello),3,'omitnan'));
iMp = squeeze(sum((Mp_B.*thkcello),3,'omitnan'));
iLp = squeeze(sum((Lp_B.*thkcello),3,'omitnan'));
iSd = squeeze(sum((Sd_B.*thkcello),3,'omitnan'));
iMd = squeeze(sum(Md_B,3,'omitnan'));
iLd = squeeze(sum(Ld_B,3,'omitnan'));
iBe = squeeze(sum(BE_B,3,'omitnan'));

%% Time series (of vert integral)
tSf = mean(iSf,1,'omitnan');
tSF = squeeze(mean(tSf,2,'omitnan'));

tSp = mean(iSp,1,'omitnan');
tSP = squeeze(mean(tSp,2,'omitnan'));

tSd = mean(iSd,1,'omitnan');
tSD = squeeze(mean(tSd,2,'omitnan'));

tMf = mean(iMf,1,'omitnan');
tMF = squeeze(mean(tMf,2,'omitnan'));

tMp = mean(iMp,1,'omitnan');
tMP = squeeze(mean(tMp,2,'omitnan'));

tLp = mean(iLp,1,'omitnan');
tLP = squeeze(mean(tLp,2,'omitnan'));

tMd = mean(iMd,1,'omitnan');
tMD = squeeze(mean(tMd,2,'omitnan'));

tLd = mean(iLd,1,'omitnan');
tLD = squeeze(mean(tLd,2,'omitnan'));

tBe = mean(iBe,1,'omitnan');
tBE = squeeze(mean(tBe,2,'omitnan'));

% tMd = mean(Md_B(:,:,35,:),1,'omitnan');
% tMd = squeeze(mean(tMd,2,'omitnan'));
% 
% tLd = mean(Ld_B(:,:,35,:),1,'omitnan');
% tLd = squeeze(mean(tLd,2,'omitnan'));
%
% tBe = mean(BE_B(:,:,35,:),1,'omitnan');
% tBE = squeeze(mean(tBe,2,'omitnan'));

%% spatial mean (of vert integral)
sSF = mean(iSf,3,'omitnan');
sSP = mean(iSp,3,'omitnan');
sSD = mean(iSd,3,'omitnan');
sMF = mean(iMf,3,'omitnan');
sMP = mean(iMp,3,'omitnan');
sLP = mean(iLp,3,'omitnan');
sMD = mean(iMd,3,'omitnan');
sLD = mean(iLd,3,'omitnan');
sBE = mean(iBe,3,'omitnan');
% sMD = mean(Md_B(:,:,35,:),4,'omitnan');
% sLD = mean(Ld_B(:,:,35,:),4,'omitnan');
% sBE = mean(BE_B(:,:,35,:),4,'omitnan');

%%
save([fpath '19900101.ocean_feisty_tracers_z_means.mat'],...
    'tSF','tSP','tSD','tMF','tMP','tMD','tLP','tLD','tBE',...
    'sSF','sSP','sSD','sMF','sMP','sMD','sLP','sLD','sBE',...
    'vSF','vSP','vSD','vMF','vMP','vLP','vMD','vLD','vBE')






