% Read NASA GISS output netcdfs

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/';

%% Each month of inputs

ncdisp([fpath '1925.Temp_bot.nc'])
%
% Source:
%            /Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/19250101.pot_temp_bot.nc
% Format:
%            classic
% Global Attributes:
%            CDI         = 'Climate Data Interface version 2.3.0 (https://mpimet.mpg.de/cdi)'
%            Conventions = 'CF-1.6'
%            xlabel      = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto      = 'From:  1925  JAN  1,  Hr  0      To:  1925  FEB  1, Hr  0  Model-Time:  2384208     Dif:  31.00 Days'
%            history     = 'Tue Oct 08 10:41:16 2024: cdo -v bottomvalue 19250101.pot_temp.nc 19250101.pot_temp_bot.nc
%                          Tue Oct 08 10:41:16 2024: cdo -selname,pot_temp 19250101.oijlVolMIPCarbANL1924Control.nc 19250101.pot_temp.nc'
%            CDO         = 'Climate Data Operators version 2.3.0 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            lono = 288
%            lato = 180
% Variables:
%     lono    
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude'
%                        units         = 'degrees_east'
%                        axis          = 'X'
%     lato    
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude'
%                        units         = 'degrees_north'
%                        axis          = 'Y'
%     pot_temp
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        long_name     = 'OCEAN POTENTIAL TEMPERATURE'
%                        units         = 'C'
%                        _FillValue    = -1.000000015047466e+30
%                        missing_value = -1.000000015047466e+30

%%
ncdisp([fpath '1925.Temp_mn_200m.nc'])

% Source:
%            /Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/19250101.pot_temp_mn_200m.nc
% Format:
%            classic
% Global Attributes:
%            CDI         = 'Climate Data Interface version 2.3.0 (https://mpimet.mpg.de/cdi)'
%            Conventions = 'CF-1.6'
%            xlabel      = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto      = 'From:  1925  JAN  1,  Hr  0      To:  1925  FEB  1, Hr  0  Model-Time:  2384208     Dif:  31.00 Days'
%            history     = 'Tue Oct 08 10:41:16 2024: cdo vertmean -select,levrange=0,200 19250101.pot_temp_bnd.nc 19250101.pot_temp_mn_200m.nc
%                          Tue Oct  8 10:41:16 2024: ncap2 -s defdim("zoc_bnd",2); zoc_bnds=make_bounds(zoc,$zoc_bnd,"zoc_bnds"); 19250101.pot_temp.nc 19250101.pot_temp_bnd.nc
%                          Tue Oct 08 10:41:16 2024: cdo -selname,pot_temp 19250101.oijlVolMIPCarbANL1924Control.nc 19250101.pot_temp.nc'
%            NCO         = 'netCDF Operators version 4.8.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)'
%            CDO         = 'Climate Data Operators version 2.3.0 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            lono = 288
%            lato = 180
% Variables:
%     lono    
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude'
%                        units         = 'degrees_east'
%                        axis          = 'X'
%     lato    
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude'
%                        units         = 'degrees_north'
%                        axis          = 'Y'
%     pot_temp
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        long_name     = 'OCEAN POTENTIAL TEMPERATURE'
%                        units         = 'C'
%                        _FillValue    = -1.000000015047466e+30
%                        missing_value = -1.000000015047466e+30

%%
ncdisp([fpath '1925.N_det_bot.nc'])

% Source:
%            /Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/19250101.N_det_bot.nc
% Format:
%            classic
% Global Attributes:
%            CDI         = 'Climate Data Interface version 2.3.0 (https://mpimet.mpg.de/cdi)'
%            Conventions = 'CF-1.6'
%            xlabel      = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto      = 'From:  1925  JAN  1,  Hr  0      To:  1925  FEB  1, Hr  0  Model-Time:  2384208     Dif:  31.00 Days'
%            history     = 'Tue Oct 08 10:41:16 2024: cdo -v bottomvalue 19250101.N_det.nc 19250101.N_det_bot.nc
%                          Tue Oct 08 10:41:16 2024: cdo -selname,N_det 19250101.toijlVolMIPCarbANL1924Control.nc 19250101.N_det.nc'
%            CDO         = 'Climate Data Operators version 2.3.0 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            lono = 288
%            lato = 180
% Variables:
%     lono 
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude'
%                        units         = 'degrees_east'
%                        axis          = 'X'
%     lato 
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude'
%                        units         = 'degrees_north'
%                        axis          = 'Y'
%     N_det
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        long_name     = 'OCEAN N_det'
%                        units         = 'kg/kg'
%                        _FillValue    = -1.000000015047466e+30
%                        missing_value = -1.000000015047466e+30

%%
ncdisp([fpath '1925.Herb_int_200m.nc'])

% Source:
%            /Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/1925.Herb_int_200m.nc
% Format:
%            classic
% Global Attributes:
%            CDI         = 'Climate Data Interface version 2.3.0 (https://mpimet.mpg.de/cdi)'
%            Conventions = 'CF-1.6'
%            xlabel      = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto      = 'From: 1925 Jan 1 hr 0  To: 1925 Dec 31 hr 0'
%            history     = 'Wed Oct 23 15:14:57 2024: cdo settaxis,1925-01-01,12:00:00,1mon 1925.Herb_int_200m.nc 1925.Herb_int_200m.nc
%                          Wed Oct 23 15:14:57 2024: ncatted -O -a fromto,global,o,c,From: 1925 Jan 1 hr 0  To: 1925 Dec 31 hr 0 1925.Herb_int_200m.nc
%                          Wed Oct 23 15:14:56 2024: ncecat -O -u time 19250101.Herb_200m.nc 19250201.Herb_200m.nc 19250301.Herb_200m.nc 19250401.Herb_200m.nc 19250501.Herb_200m.nc 19250601.Herb_200m.nc 19250701.Herb_200m.nc 19250801.Herb_200m.nc 19250901.Herb_200m.nc 19251001.Herb_200m.nc 19251101.Herb_200m.nc 19251201.Herb_200m.nc 1925.Herb_int_200m.nc
%                          Wed Oct 23 15:14:36 2024: cdo setattribute,Herb_wt@long_name=HERB INTEGRATED OVER TOP 200m 19250101.Herb_200m.nc 19250101.Herb_200m.nc
%                          Wed Oct 23 15:14:36 2024: cdo setattribute,Herb_wt@units=kg/m2 19250101.Herb_200m.nc 19250101.Herb_200m.nc
%                          Wed Oct 23 15:14:36 2024: cdo vertsum -select,levrange=0,200 tempfile2.nc 19250101.Herb_200m.nc
%                          Wed Oct 23 15:14:36 2024: cdo expr,Herb_wt=mo*Herb tempfile1.nc tempfile2.nc
%                          Wed Oct 23 15:14:36 2024: cdo merge 19250101.VolFluid.nc 19250101.Herb_bnd.nc tempfile1.nc
%                          Wed Oct 23 15:14:34 2024: cdo -selname,mo 19250101.oijlVolMIPCarbANL1924Control.nc 19250101.VolFluid.nc'
%            NCO         = 'netCDF Operators version 4.8.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)'
%            CDO         = 'Climate Data Operators version 2.3.0 (https://mpimet.mpg.de/cdo)'
% Dimensions:
%            time = 12    (UNLIMITED)
%            lono = 288
%            lato = 180
% Variables:
%     time   
%            Size:       12x1
%            Dimensions: time
%            Datatype:   double
%            Attributes:
%                        standard_name = 'time'
%                        units         = 'month as %Y%m.%f'
%                        calendar      = 'proleptic_gregorian'
%                        axis          = 'T'
%     lono   
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        standard_name = 'longitude'
%                        long_name     = 'longitude'
%                        units         = 'degrees_east'
%                        axis          = 'X'
%     lato   
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        standard_name = 'latitude'
%                        long_name     = 'latitude'
%                        units         = 'degrees_north'
%                        axis          = 'Y'
%     Herb_wt
%            Size:       288x180x12
%            Dimensions: lono,lato,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'HERB INTEGRATED OVER TOP 200m'
%                        units         = 'kg/m2'
%                        _FillValue    = -8.999999873090293e+33
%                        missing_value = -8.999999873090293e+33

%%
ncid = netcdf.open([fpath '19250101.Herb_mn_100m.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

Herb100 = Herb;

%% 150
clear Herb

ncid = netcdf.open([fpath '19250101.Herb_mn_150m.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

Herb150 = Herb;

%% 200
clear Herb

ncid = netcdf.open([fpath '19250101.Herb_mn_200m.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

Herb200 = Herb;

%%
clear Herb

%%
Herb100 = double(Herb100);
Herb150 = double(Herb150);
Herb200 = double(Herb200);

Herb100(Herb100<= -1.00e+30) = nan;
Herb150(Herb150<= -1.00e+30) = nan;
Herb200(Herb200<= -1.00e+30) = nan;

%%
figure(1)
pcolor(log10(Herb100)); shading flat; colorbar;
clim([-12 -7])

figure(2)
pcolor(log10(Herb150)); shading flat; colorbar;
clim([-12 -7])

figure(3)
pcolor(log10(Herb200)); shading flat; colorbar;
clim([-12 -7])

%%
%save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat']);




