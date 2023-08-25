% Read NASA GISS output netcdfs

clear 
close all

fpath='/Volumes/petrik-lab/Feisty/GCM_Data/GISS/VolMIP/';

%% Each month of inputs
ncdisp([fpath 'APR1925.aijVolMIPCarbANL1924Control.nc'])
% all atmo

%%
ncdisp([fpath 'APR1925.obio_ijlVolMIPCarbANL1924Control.nc'])
% ocean bio 3D vars (i,j,l = x,y,z)

% Global Attributes:
%            xlabel = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto = 'From:  1925  APR  1,  Hr  0      To:  1925  MAY  1, Hr  0  Model-Time:  2388480     Dif:  30.00 Days'
% Dimensions:
%            lono  = 288
%            lono2 = 288
%            lato  = 180
%            lato2 = 180
%            zoc   = 40
%            zoce  = 40
% Variables:
%     lono     
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_east'
%     lono2    
%            Size:       288x1
%            Dimensions: lono2
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_east'
%     lato     
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_north'
%     lato2    
%            Size:       180x1
%            Dimensions: lato2
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_north'
%     zoc      
%            Size:       40x1
%            Dimensions: zoc
%            Datatype:   single
%            Attributes:
%                        units    = 'm'
%                        positive = 'down'
%     zoce     
%            Size:       40x1
%            Dimensions: zoce
%            Datatype:   single
%            Attributes:
%                        units    = 'm'
%                        positive = 'down'
%     avgq     
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'quanta/m2/s'
%                        long_name     = 'Mean daily irradiance'
%                        missing_value = -1.000000015047466e+30
%     diatwss  
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'diatwss'
%                        missing_value = -1.000000015047466e+30
%     chlowss  
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'chlowss'
%                        missing_value = -1.000000015047466e+30
%     cyanwss  
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'cyanwss'
%                        missing_value = -1.000000015047466e+30
%     coccwss  
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'coccwss'
%                        missing_value = -1.000000015047466e+30
%     ndetwsdet
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'ndetwsdet'
%                        missing_value = -1.000000015047466e+30
%     sdetwsdet
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'sdetwsdet'
%                        missing_value = -1.000000015047466e+30
%     idetwsdet
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'm/s'
%                        long_name     = 'idetwsdet'
%                        missing_value = -1.000000015047466e+30
%     diatpp   
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'mg,C/m2/day'
%                        long_name     = 'diatpp'
%                        missing_value = -1.000000015047466e+30
%     chlopp   
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'mg,C/m2/day'
%                        long_name     = 'chlopp'
%                        missing_value = -1.000000015047466e+30
%     cyanpp   
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'mg,C/m2/day'
%                        long_name     = 'cyanpp'
%                        missing_value = -1.000000015047466e+30
%     coccpp   
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'mg,C/m2/day'
%                        long_name     = 'coccpp'
%                        missing_value = -1.000000015047466e+30
%     diat_lim1,2,3,4,5
%     chlo_lim1,2,3,4,5
%     cyan_lim1,2,3,4,5
%     cocc_lim1,2,3,4,5
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = '?'
%                        long_name     = 'limiting factor #'
%                        missing_value = -1.000000015047466e+30

%%
ncdisp([fpath 'APR1925.obio_ijVolMIPCarbANL1924Control.nc'])
% ocean BGC 2D vars (i,j = x,y)

% Global Attributes:
%            xlabel                    = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto                    = 'From:  1925  APR  1,  Hr  0      To:  1925  MAY  1, Hr  0  Model-Time:  2388480     Dif:  30.00 Days'
%            history                   = 'Fri May 27 10:11:49 2022: ncks -A -v oxyp /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc APR1925.obio_ijVolMIPCarbANL1924Control.nc
%                                        Fri May 27 10:07:36 2022: ncks -A -v oxyp /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc APR1925.obio_ijVolMIPCarbANL1924Control.nc
%                                        Tue Mar 22 14:52:25 2022: ncks -A -v oxyp /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc APR1925.obio_ijVolMIPCarbANL1924Control.nc
%                                        Mon Mar 21 13:12:29 2022: ncks -A -v oxyp /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc APR1925.obio_ijVolMIPCarbANL1924Control.nc
%                                        Mon Mar 21 00:13:39 2022: ncks -A -v oxyp /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc APR1925.obio_ijVolMIPCarbANL1924Control.nc'
%            history_of_appended_files = 'Fri May 27 10:11:49 2022: Appended file /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc had following "history" attribute:
%                                        Tue Sep 14 15:08:01 2021: ncks -v oxyp SEP1973.oijVolMIPCarbANL1968c.nc VolMIPCarbOXYP.nc
%                                        Fri May 27 10:07:36 2022: Appended file /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc had following "history" attribute:
%                                        Tue Sep 14 15:08:01 2021: ncks -v oxyp SEP1973.oijVolMIPCarbANL1968c.nc VolMIPCarbOXYP.nc
%                                        Tue Mar 22 14:52:25 2022: Appended file /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc had following "history" attribute:
%                                        Tue Sep 14 15:08:01 2021: ncks -v oxyp SEP1973.oijVolMIPCarbANL1968c.nc VolMIPCarbOXYP.nc
%                                        Mon Mar 21 13:12:29 2022: Appended file /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc had following "history" attribute:
%                                        Tue Sep 14 15:08:01 2021: ncks -v oxyp SEP1973.oijVolMIPCarbANL1968c.nc VolMIPCarbOXYP.nc
%                                        Mon Mar 21 00:13:39 2022: Appended file /archive/u/dlitchmo/VolmipProjectHub/RawScaled/suppData/VolMIPCarbOXYP.nc had following "history" attribute:
%                                        Tue Sep 14 15:08:01 2021: ncks -v oxyp SEP1973.oijVolMIPCarbANL1968c.nc VolMIPCarbOXYP.nc
%                                        '
%            NCO                       = 'netCDF Operators version 4.8.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)'
% Dimensions:
%            lono  = 288
%            lono2 = 288
%            lato  = 180
%            lato2 = 180
% Variables:
%     lono         
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_east'
%     lono2        
%            Size:       288x1
%            Dimensions: lono2
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_east'
%     lato         
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_north'
%     lato2        
%            Size:       180x1
%            Dimensions: lato2
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_north'
%     oij_pH       
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'pH units'
%                        long_name = 'ocean surface pH'
%     oij_dayl     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'timesteps'
%                        long_name = 'Daylight length'
%     oij_diat     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg/m3'
%                        long_name = 'Surface ocean Diatoms'
%     oij_chlo     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg/m3'
%                        long_name = 'Surface ocean Chlorophytes'
%     oij_cyan     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg/m3'
%                        long_name = 'Surface ocean Cyanobacteria'
%     oij_cocc     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg/m3'
%                        long_name = 'Surface ocean Coccolithophores'
%     oij_herb     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg/m3'
%                        long_name = 'Surface ocean Herbivores'
%     oij_flux     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'depends if on atm/ocean grid'
%                        long_name = 'AO Flux CO2 (gr,CO2 or mol,CO2/m2/yr)'
%     oij_cexp     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'PgC/yr'
%                        long_name = 'C export flux at compensation depth'
%     oij_ndet     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'ugC/l'
%                        long_name = 'N/C detritus at 74m'
%     oij_setl     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'm/s'
%                        long_name = 'settlvel n/cdet at 74m'
%     oij_sink     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'm/s'
%                        long_name = 'sink vel phytopl at 74m'
%     oij_xchl     
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'kg,C*m/s'
%                        long_name = 'C export due to chloroph'
%     oij_fca      
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mili-g,C/m2/s'
%                        long_name = 'CaCO3 export flux at compensation depth'
%     oij_pp       
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg,C/m2/day'
%                        long_name = 'Depth integrated PP'
%     oij_pp1      
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg,C/m2/day'
%                        long_name = 'PP-diat'
%     oij_pp2      
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg,C/m2/day'
%                        long_name = 'PP-chlor'
%     oij_pp3      
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg,C/m2/day'
%                        long_name = 'PP-cyan'
%     oij_pp4      
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = 'mg,C/m2/day'
%                        long_name = 'PP-cocc'

%%
ncdisp([fpath 'APR1925.oijlVolMIPCarbANL1924Control.nc'])
% ocean 3D velocities and fluxes

% pot_temp    
%            Size:       288x180x40
%            Dimensions: lono,lato,zoc
%            Datatype:   single
%            Attributes:
%                        units         = 'C'
%                        long_name     = 'OCEAN POTENTIAL TEMPERATURE'
%                        missing_value = -1.000000015047466e+30
%     mfw2        
%            Size:       288x180x40
%            Dimensions: lono,lato,zoce
%            Datatype:   single
%            Attributes:
%                        units         = 'kg^2/m^4'
%                        long_name     = 'Ocean vertical mass flux squared'
%                        missing_value = -1.000000015047466e+30

%%
ncdisp([fpath 'APR1925.oijVolMIPCarbANL1924Control.nc'])
%ocean 2D variables

% Global Attributes:
%            xlabel = 'VolMIPCarbANL1924Control (VolMIPCarbANL1924 + updated aerosol/ozone input files'
%            fromto = 'From:  1925  APR  1,  Hr  0      To:  1925  MAY  1, Hr  0  Model-Time:  2388480     Dif:  30.00 Days'
% Dimensions:
%            lono  = 288
%            lono2 = 288
%            lato  = 180
%            lato2 = 180
% Variables:
%     lono       
%            Size:       288x1
%            Dimensions: lono
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_east'
%     lono2      
%            Size:       288x1
%            Dimensions: lono2
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_east'
%     lato       
%            Size:       180x1
%            Dimensions: lato
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_north'
%     lato2      
%            Size:       180x1
%            Dimensions: lato2
%            Datatype:   single
%            Attributes:
%                        units = 'degrees_north'
%     oij_mask   
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units     = '1'
%                        long_name = 'Ocean Mask'
%     oij_mld    
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units         = 'm'
%                        long_name     = 'Ocean Mixed layer depth'
%                        missing_value = -1.000000015047466e+30
%     oij_ssh    
%            Size:       288x180
%            Dimensions: lono,lato
%            Datatype:   single
%            Attributes:
%                        units         = 'm'
%                        long_name     = 'Ocean surface height'
%                        missing_value = -1.000000015047466e+30

%%
ncdisp([fpath 'APR1925.ojlVolMIPCarbANL1924Control.nc'])
% ocean basin-scale variables

%%
ncdisp([fpath 'APR1925.toijlVolMIPCarbANL1924Control.nc'])
% 3D ocean tracers (all fluxes)

%%
ncid = netcdf.open([fpath 'APR1925.oijlVolMIPCarbANL1924Control.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:(nvars)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1.000000020040877e+20) = NaN;']);
end
netcdf.close(ncid);

%% Time
yr = (time-time(1)+1)/365;

%%
save([fpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat']);





