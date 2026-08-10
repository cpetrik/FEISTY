% 0.5 degree global MOM6-COBALTv3 
% min vol of grid cell


clear
close all

fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/OM4_05_COBALTv3_FEISTYoff/';

%%
ni=720;
nj=576;
nk=35;
nt=60;

%% vertical
ncid = netcdf.open([fpath 'ocean_cobalt_feisty_forcing_z.199001-199412.volcello.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 8
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1, [0,0,0,0],[ni nj nk 12]);']);
    eval([ varname '(' varname ' == 1e20) = NaN;']);
end

%%
volcello(volcello>1e19)=nan;
mvol = min(volcello(:,:,1:10,:),[],4,'omitnan');

%%
figure()
subplot(2,2,1)
pcolor(squeeze(mvol(:,:,1))); shading flat; colorbar
subplot(2,2,2)
pcolor(squeeze(mvol(:,:,2))); shading flat; colorbar
subplot(2,2,3)
pcolor(squeeze(mvol(:,:,3))); shading flat; colorbar
subplot(2,2,4)
pcolor(squeeze(mvol(:,:,4))); shading flat; colorbar

%%
figure()
histogram(mvol(:))

mean(mvol(:),'omitnan') % 3.0183e+10 m3

%%
zmvol = squeeze(mean(mvol,[1 2],'omitnan'));

%%
0.02/3e10

11.2/3e10

5.6e3/3e10