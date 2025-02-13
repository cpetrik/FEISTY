% Read NASA GISS output netcdfs

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/';

%%
ncdisp([fpath 'Data2Colleen.nc'])

%%
ncid = netcdf.open([fpath 'Data2Colleen.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Rearrange
% test = ones(3,3,12,10);
% for m=1:12
%     for y=1:10
%     test(:,:,m,y) = m*y;
%     end
% end
% test2=reshape(test,3,3,12*10,1);

[ni,nj,nm,ny] = size(ptemp_upper200m_avg);
ptemp_upper200m_avg = reshape(ptemp_upper200m_avg,ni,nj,nm*ny);

% mpt=mean(ptemp_upper200m_avg,1,'omitnan');
% mpt=mean(mpt,2,'omitnan');
% mpt=squeeze(mpt);
% plot(mpt)

ptemp_bottom = reshape(ptemp_bottom,ni,nj,nm*ny);
poc_mgcm3_bottom = reshape(poc_mgcm3_bottom,ni,nj,nm*ny);
biomass_Chlo_intrg200m = reshape(biomass_Chlo_intrg200m,ni,nj,nm*ny);
biomass_Cocc_intrg200m = reshape(biomass_Cocc_intrg200m,ni,nj,nm*ny);
biomass_Cyan_intrg200m = reshape(biomass_Cyan_intrg200m,ni,nj,nm*ny);
biomass_Diat_intrg200m = reshape(biomass_Diat_intrg200m,ni,nj,nm*ny);
biomass_Herb_intrg200m = reshape(biomass_Herb_intrg200m,ni,nj,nm*ny);

%% L frac
Lfrac = biomass_Diat_intrg200m ./ (biomass_Diat_intrg200m+biomass_Chlo_intrg200m+...
    biomass_Cocc_intrg200m+biomass_Cyan_intrg200m);
biomass_Zmeso_intrg200m = Lfrac .* biomass_Herb_intrg200m;

%%
save([fpath 'giss_10yr_FEISTY_forcing.mat'],'biomass_Diat_intrg200m',...
    'biomass_Chlo_intrg200m','biomass_Cocc_intrg200m','biomass_Cyan_intrg200m',...
    'biomass_Herb_intrg200m','biomass_Zmeso_intrg200m','Lfrac','poc_mgcm3_bottom',...
    'ptemp_bottom','ptemp_upper200m_avg');

save([fpath 'giss_10yr_gridspec.mat'],'lat','lon');
