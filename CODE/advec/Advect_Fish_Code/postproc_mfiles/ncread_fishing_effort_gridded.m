% Jerome's GFW data. Gridding of the last update of the fishing effort
% in h per 1deg grid cell, monthly, over 10 years (2015-2024).
% I have cut out the first 2 years because data coverage had a lot of gaps
% There are 16 different gears (see names in the .nc metadata),
% usually I use [1,4,8,14] for demersal fishing
% [2,5,6,7,9,10,13,16] for pelagic fishing.
% Note that category (3) 'fishing' is, if I recall properly, a category
% that includes fishing activity that was not clearly labelled.
% I usually disregard it for the pel vs. dem analysis.

clear
close all

%%
epath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/movement/FishingDataColleen/';

%%
ncdisp([epath 'Global_fishing_effort_1x1_per_gear_2015-2024.nc'])

%Dimensions:
%            lon   = 360
%            lat   = 180
%            year  = 10    (UNLIMITED)
%            month = 12
%            gear  = 16
% Variables:
%     lon
%            Size:       360x1
%            Dimensions: lon
%            Datatype:   double
%            Attributes:
%                        long_name = 'longitude'
%                        units     = 'degrees_east'
%                        axis      = 'X'
%     lat
%            Size:       180x1
%            Dimensions: lat
%            Datatype:   double
%            Attributes:
%                        long_name = 'latitude'
%                        units     = 'degrees_north'
%                        axis      = 'Y'
%     year
%            Size:       10x1
%            Dimensions: year
%            Datatype:   double
%            Attributes:
%                        long_name = 'year'
%     month
%            Size:       12x1
%            Dimensions: month
%            Datatype:   double
%            Attributes:
%                        long_name = 'month'
%                        axis      = 'M'
%     gear
%            Size:       16x1
%            Dimensions: gear
%            Datatype:   double
%            Attributes:
%                        long_name = '(1) dredge_fishing; (2) drifting_longlines; (3) fishing; (4) fixed_gear; (5) other_purse_seines; (6) other_seines; (7) pole_and_line; (8) pots_and_traps; (9) purse_seines; (10) seiners; (11) set_gillnets; (12) set_longlines; (13) squid_jigger; (14) trawlers; (15) trollers; (16) tuna_purse_seines'
%     effort_h
%            Size:       360x180x12x16x10
%            Dimensions: lon,lat,month,gear,year
%            Datatype:   single
%            Attributes:
%                        _FillValue    = NaN
%                        missing_value = 1e+20
%                        units         = 'h'

%%
gear_text = {'dredge_fishing','drifting_longlines','fishing','fixed_gear',...
    'other_purse_seines','other_seines','pole_and_line','pots_and_traps',...
    'purse_seines','seiners','set_gillnets','set_longlines','squid_jigger',...
    'trawlers','trollers','tuna_purse_seines'}';

%% 
ncid = netcdf.open([epath 'Global_fishing_effort_1x1_per_gear_2015-2024.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 1e+20) = NaN;']);
end

%% get rid of nans 
effort_h(effort_h>1e19) = nan;

%% gears
%[1,4,8,14] for demersal fishing
% [2,5,6,7,9,10,13,16] for pelagic fishing.
did = [1,4,8,14];
pid = [2,5,6,7,9:13,16];

Deffort = effort_h(:,:,:,did,:);
Peffort = effort_h(:,:,:,pid,:);

sumDeffort = squeeze(sum(Deffort,4,'omitnan'));
sumPeffort = squeeze(sum(Peffort,4,'omitnan'));

sumDeffort = reshape(sumDeffort,360,180,120);
sumPeffort = reshape(sumPeffort,360,180,120);

%% means over full time
meanDeffort = double(squeeze(mean(sumDeffort,3,'omitnan')));
meanPeffort = double(squeeze(mean(sumPeffort,3,'omitnan')));

%%
save([epath 'Global_fishing_effort_1deg_pel_dem_mean2015-2024.mat'],...
    'meanDeffort','meanPeffort','did','pid','gear_text','lat','lon');


