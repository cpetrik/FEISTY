%Example of dividing zooplankton into multiple size bins
%Daewel et al 2008 method

%Mesozoop biom ranges from 1e2 to 1e4 mg C/m2

%BMZ = 3.6; 5.8 Pg C/yr =1e15gC/yr
%BLZ = 1.1; 1.7 Pg C/yr

%%
%365d/yr
%360e6 km2

BMZ = 3.6e15 / 365 / 360e6; %/ 1e6; %in gC/km2
BLZ = 1.1e15 / 365 / 360e6; %/ 1e6;

%%

mlow = 0.2;
mhi = 2;
llow = 2; 
lhi = 20;
hm = 0.2;
hl = 2;
mbin_ends = mlow:hm:mhi;
mbin_mids = ((mbin_ends(1)+mbin_ends(2))/2):hm:((mbin_ends(end-1)+mbin_ends(end))/2);
lbin_ends = llow:hl:lhi;
lbin_mids = ((lbin_ends(1)+lbin_ends(2))/2):hl:((lbin_ends(end-1)+lbin_ends(end))/2);
Hm = mbin_ends(end) - mbin_ends(1);
Hl = lbin_ends(end) - lbin_ends(1);
Mmed = (mlow+mhi)/2;
Lmed = (llow+lhi)/2;

%% size in um
pl = 1e3* [mbin_mids,lbin_mids]';

%% size distr
SD = 695.73 .* exp(-0.0083.*pl);

%% ind mass in ug C
mug = exp(2.772 * log(pl) - 7.476);
mg = 1e-6*mug;

%% total biomass in each size bin
tm = sum(SD.*mg);
SP = SD * mg / tm;

%% solve for SD given tm and fixed pl & m
n=length(pl);
tm=BMZ+BLZ;


