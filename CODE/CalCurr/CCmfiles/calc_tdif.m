% Calculate fraction pelagic vs. benthic for dems

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
%figp = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
sname = 'Clim_';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);
load([dpath sname 'locs_' harv '.mat'])

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%%
SP = S_Sml_p;
SF = S_Sml_f;
SD = S_Sml_d;
MP = S_Med_p;
MF = S_Med_f;
MD = S_Med_d;
LP = S_Lrg_p;
LD = S_Lrg_d;
CO = S_Cobalt;

t=1:size(SP,1);
lyr=t((end-12+1):end);

%% Final mean biomass in each size
Sm = 0.25;  %Feeding 2 sizes down
J = 1.0;    %Juvenile feeding reduction
D = 0.75;   %Demersal feeding in pelagic reduction
A = 0.50;

LD_phi_MF = D*A;
LD_phi_MP = D;
LD_phi_MD = 1.0;
LD_phi_BE = 1.0;

bio1=squeeze((MF(:,1,:)));
bio2=squeeze((MP(:,1,:)));
bio3=squeeze((MD(:,1,:)));
bio4=squeeze((CO(:,1,:)));

biop = LD_phi_MF*bio1 + LD_phi_MP*bio2;
biod = LD_phi_MD*bio3 + LD_phi_BE*bio4;

tdif = biop ./ (biop + biod);

%% means in last year
coast = [3;4;6;8;10];
coastal=spots(coast)';
mtdiff=mean(tdif(lyr,coast));
save([dpath sname 'Coastal_locs_LD_frac_time_pel.png'],'tdif','mtdiff',...
    'spots','coast','coastal');

%%
yr=t/12;
figure(1)
plot(yr,tdif(:,coast),'LineWidth',1.5)
legend(spots{coast})
%legend('location','eastoutside')
xlabel('time (yrs)')
ylabel('fraction time pelagic')
stamp(cfile)
print('-dpng',[fpath sname  'Coastal_LD_frac_time_pel.png'])

