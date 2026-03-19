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
%ncdisp([fpath '19900101.ocean_feisty_pelagic_fluxes_z.nc'])

%% want met, prod, E_A, f_tot, Fout, rho


st = 1990;

tLD_met2 = nan*ones(1,12);
tLD_prod2 = tLD_met2;
tLD_rho2 = tLD_met2;
tLP_met2 = tLD_met2;
tLP_prod2 = tLD_met2;
tLP_rho2 = tLD_met2;
tMD_met2 = tLD_met2;
tMD_Fout2 = tLD_met2;
tMD_prod2 = tLD_met2;
tMP_met2 = tLD_met2;
tMP_Fout2 = tLD_met2;
tMP_prod2 = tLD_met2;
tMF_met2 = tLD_met2;
tMF_prod2 = tLD_met2;
tMF_rho2 = tLD_met2;
tSD_met2 = tLD_met2;
tSD_Fout2 = tLD_met2;
tSD_prod2 = tLD_met2;
tSP_met2 = tLD_met2;
tSP_Fout2 = tLD_met2;
tSP_prod2 = tLD_met2;
tSF_met2 = tLD_met2;
tSF_Fout2 = tLD_met2;
tSF_prod2 = tLD_met2;
tMD_EA2 = tLD_met2;
tMP_EA2 = tLD_met2;
tMF_EA2 = tLD_met2;
tSD_EA2 = tLD_met2;
tSP_EA2 = tLD_met2;
tSF_EA2 = tLD_met2;

mLD_met2 = nan*ones(720,576,length(st));
mLD_prod2 = mLD_met2;
mLD_rho2 = mLD_met2;
mLP_met2 = mLD_met2;
mLP_prod2 = mLD_met2;
mLP_rho2 = mLD_met2;
mMD_met2 = mLD_met2;
mMD_Fout2 = mLD_met2;
mMD_prod2 = mLD_met2;
mMP_met2 = mLD_met2;
mMP_Fout2 = mLD_met2;
mMP_prod2 = mLD_met2;
mMF_met2 = mLD_met2;
mMF_prod2 = mLD_met2;
mMF_rho2 = mLD_met2;
mSD_met2 = mLD_met2;
mSD_Fout2 = mLD_met2;
mSD_prod2 = mLD_met2;
mSP_met2 = mLD_met2;
mSP_Fout2 = mLD_met2;
mSP_prod2 = mLD_met2;
mSF_met2 = mLD_met2;
mSF_Fout2 = mLD_met2;
mSF_prod2 = mLD_met2;
mMD_EA2 = mLD_met2;
mMP_EA2 = mLD_met2;
mMF_EA2 = mLD_met2;
mSD_EA2 = mLD_met2;
mSP_EA2 = mLD_met2;
mSF_EA2 = mLD_met2;

vLD_met2 = nan*ones(35,length(st));
vLD_prod2 = vLD_met2;
vLD_rho2 = vLD_met2;
vLP_met2 = vLD_met2;
vLP_prod2 = vLD_met2;
vLP_rho2 = vLD_met2;
vMD_met2 = vLD_met2;
vMD_Fout2 = vLD_met2;
vMD_prod2 = vLD_met2;
vMP_met2 = vLD_met2;
vMP_Fout2 = vLD_met2;
vMP_prod2 = vLD_met2;
vMF_met2 = vLD_met2;
vMF_prod2 = vLD_met2;
vMF_rho2 = vLD_met2;
vSD_met2 = vLD_met2;
vSD_Fout2 = vLD_met2;
vSD_prod2 = vLD_met2;
vSP_met2 = vLD_met2;
vSP_Fout2 = vLD_met2;
vSP_prod2 = vLD_met2;
vSF_met2 = vLD_met2;
vSF_Fout2 = vLD_met2;
vSF_prod2 = vLD_met2;
vMD_EA2 = vLD_met2;
vMP_EA2 = vLD_met2;
vMF_EA2 = vLD_met2;
vSD_EA2 = vLD_met2;
vSP_EA2 = vLD_met2;
vSF_EA2 = vLD_met2;

tMD_flev = tLD_met2;
tMP_flev = tLD_met2;
tMF_flev = tLD_met2;
tSD_flev = tLD_met2;
tSP_flev = tLD_met2;
tSF_flev = tLD_met2;
mMD_flev = mLD_met2;
mMP_flev = mLD_met2;
mMF_flev = mLD_met2;
mSD_flev = mLD_met2;
mSP_flev = mLD_met2;
mSF_flev = mLD_met2;
vMD_flev = vLD_met2;
vMP_flev = vLD_met2;
vMF_flev = vLD_met2;
vSD_flev = vLD_met2;
vSP_flev = vLD_met2;
vSF_flev = vLD_met2;

%%
ncid = netcdf.open([fpath '19900101.ocean_feisty_tracers_z.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
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
BE_B(BE_B>1e19) = nan;

%% thkcello
load([gpath 'ocean_cobalt_feisty_forcing_z.199001-199412.thkcello.mat'])

%% Vertical means
vSf = mean(Sf_B,1,'omitnan');
vSf = squeeze(mean(vSf,2,'omitnan'));
vSf = squeeze(mean(vSf,2,'omitnan'));

vSp = mean(Sp_B,1,'omitnan');
vSp = squeeze(mean(vSp,2,'omitnan'));
vSp = squeeze(mean(vSp,2,'omitnan'));

vSd = mean(Sd_B,1,'omitnan');
vSd = squeeze(mean(vSd,2,'omitnan'));
vSd = squeeze(mean(vSd,2,'omitnan'));

vMf = mean(Mf_B,1,'omitnan');
vMf = squeeze(mean(vMf,2,'omitnan'));
vMf = squeeze(mean(vMf,2,'omitnan'));

vMp = mean(Mp_B,1,'omitnan');
vMp = squeeze(mean(vMp,2,'omitnan'));
vMp = squeeze(mean(vMp,2,'omitnan'));

vLp = mean(Lp_B,1,'omitnan');
vLp = squeeze(mean(vLp,2,'omitnan'));
vLp = squeeze(mean(vLp,2,'omitnan'));

vLd = mean(Ld_B,1,'omitnan');
vLd = squeeze(mean(vLd,2,'omitnan'));
vLd = squeeze(sum(vLd,2,'omitnan'));


%% vert sums
iSf = squeeze(sum((Sf_B.*thkcello),3,'omitnan'));
iSp = squeeze(sum((Sp_B.*thkcello),3,'omitnan'));
iSd = squeeze(sum((Sd_B.*thkcello),3,'omitnan'));
iMf = squeeze(sum((Mf_B.*thkcello),3,'omitnan'));
iMp = squeeze(sum((Mp_B.*thkcello),3,'omitnan'));
iLp = squeeze(sum((Lp_B.*thkcello),3,'omitnan'));

%% Time series (of vert integral)
tSf = mean(iSf,1,'omitnan');
tSf = squeeze(mean(tSf,2,'omitnan'));

tSp = mean(iSp,1,'omitnan');
tSp = squeeze(mean(tSp,2,'omitnan'));

tSd = mean(iSd,1,'omitnan');
tSd = squeeze(mean(tSd,2,'omitnan'));

tMf = mean(iMf,1,'omitnan');
tMf = squeeze(mean(tMf,2,'omitnan'));

tMp = mean(iMp,1,'omitnan');
tMp = squeeze(mean(tMp,2,'omitnan'));

tLp = mean(iLp,1,'omitnan');
tLp = squeeze(mean(tLp,2,'omitnan'));

tMd = mean(Md_B(:,:,35,:),1,'omitnan');
tMd = squeeze(mean(tMd,2,'omitnan'));

tLd = mean(Ld_B(:,:,35,:),1,'omitnan');
tLd = squeeze(mean(tLd,2,'omitnan'));

tBe = mean(BE_B(:,:,35,:),1,'omitnan');
tBe = squeeze(mean(tBe,2,'omitnan'));

%% spatial mean (of vert integral)
sSf = mean(iSf,3,'omitnan');
sSp = mean(iSp,3,'omitnan');
sSd = mean(iSd,3,'omitnan');
sMf = mean(iMf,3,'omitnan');
sMp = mean(iMp,3,'omitnan');
sLp = mean(iLp,3,'omitnan');

sMd = mean(Md_B(:,:,35,:),4,'omitnan');
sLd = mean(Ld_B(:,:,35,:),4,'omitnan');
sBe = mean(BE_B(:,:,35,:),4,'omitnan');

%% put in arrays
y=1;
yid = (((y-1)*60)+1):(y*60);

tSF(1,yid) = tSf;
tSP(1,yid) = tSp;
tSD(1,yid) = tSd;
tMF(1,yid) = tMf;
tMP(1,yid) = tMp;
tLP(1,yid) = tLp;
tMD(1,yid) = tMd';
tLD(1,yid) = tLd';
tLP(1,yid) = tBe';

sSF(:,:,y) = sSf;
sSP(:,:,y) = sSp;
sSD(:,:,y) = sSd;
sMF(:,:,y) = sMf;
sMP(:,:,y) = sMp;
sMD(:,:,y) = sMd;
sLP(:,:,y) = sLp;
sLD(:,:,y) = sLd;
sBE(:,:,y) = sBe;

vSF(:,y) = vSf;
vSP(:,y) = vSp;
vSD(:,y) = vSd;
vMF(:,y) = vMf;
vMP(:,y) = vMp;
vLP(:,y) = vLp;

save([fpath '19900101.ocean_feisty_fluxes_z_means.mat'],'tSF','tSP','tSD',...
    'tMF','tMP','tMD','tLP','tLD','tBE','sSF','sSP','sSD',...
    'sMF','sMP','sMD','sLP','sLD','sBE','vSF','vSP','vSD',...
    'vMF','vMP','vLP')






