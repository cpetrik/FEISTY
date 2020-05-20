% Visualize production of F and P as fn of zoop
% at 2 different temps (or bins)
% ESM2M Hindcast of 1951-2000 and Forecast of 2051-2100

clear all
close all

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
global bent_eff rfrac D J Sm A benc bcmx amet
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel Jsel efn cfn mfn
global dfrate frate

cfn=nan;
efn=nan;
mfn=nan;

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
bpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% Hindcast grid
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t'); %geolon_t,geolat_t
grid = csvread([cpath 'grid_csv.csv']); %grid
ID = grid(:,1);

%% Zoop biomass
load([bpath 'cobalt_det_temp_zoop_npp_means.mat'],'mzprod_mean_hist',...
    'lzprod_mean_hist','mzprod_mean_fore','lzprod_mean_fore','ptemp_mean_hist',...
    'ptemp_mean_fore');

% molN/m2/s --> g/m2/d
mzprod_hist = mzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_hist = lzprod_mean_hist * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
mzprod_fore = mzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;
lzprod_fore = lzprod_mean_fore * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 * 24;

z_hist = mzprod_hist + lzprod_hist;
z_fore = mzprod_fore + lzprod_fore;

vhT = ptemp_mean_hist(ID);
vhZ = z_hist(ID);
vcT = ptemp_mean_fore(ID);
vcZ = z_fore(ID);

% mean temps around 10 and 25 C
hlid = find(vhT<11 & vhT>=9);
clid = find(vcT<11 & vcT>=9);

hhid = find(vhT<26 & vhT>=24);
chid = find(vcT<26 & vcT>=24);


%% FEISTY Output
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

harv = 'All_fish03';

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
load([nfile 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params');
params = red_params;

%%
epred = nan*ones(13,4,length(params));
for j = 1:length(params)
    %! Change individual parameters
    pset = params(j,:);
    set_params6_samek(pset)
    
    %! Make core parameters/constants (global)
    const_params6_samek()
    
    %% Historic
    %! Create a directory for output
    [fname,simname] = sub_fname_hist_ens(frate);
        
    % Last 50 year means
    load([fname '_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50',...
        'mf_prod50','mp_prod50',...
        'lp_prod50');
    hF = sf_prod50+mf_prod50;
    hP = sp_prod50+mp_prod50+lp_prod50;
    
    %% Forecast
    %! Create a directory for output
    [fname,simname] = sub_fname_fore_ens(frate);
        
    % Last 50 year means
    load([fname '_Means_prod_' simname '.mat'],...
        'sf_prod50','sp_prod50',...
        'mf_prod50','mp_prod50',...
        'lp_prod50');
    cF = sf_prod50+mf_prod50;
    cP = sp_prod50+mp_prod50+lp_prod50;



%%

vhF = hF;
vhP = hP;

vcF = cF;
vcP = cP;

%%
% figure(1)
% subplot(2,2,1)
% plot(vhZ(hlid),vhP(hlid),'bo'); hold on;
% plot(vhZ(hlid),vhF(hlid),'ro'); hold on;
% title('10C hist')
% axis([0 25 0 20])
% 
% subplot(2,2,2)
% plot(vhZ(clid),vhP(clid),'bo'); hold on;
% plot(vhZ(clid),vhF(clid),'ro'); hold on;
% title('10C fore')
% axis([0 25 0 20])
% 
% subplot(2,2,3)
% plot(vhZ(hhid),vhP(hhid),'bo'); hold on;
% plot(vhZ(hhid),vhF(hhid),'ro'); hold on;
% title('25C hist')
% axis([0 25 0 20])
% 
% subplot(2,2,4)
% plot(vhZ(chid),vhP(chid),'bo'); hold on;
% plot(vhZ(chid),vhF(chid),'ro'); hold on;
% title('25C fore')
% axis([0 25 0 20])

%% combine past anf future, fit lines
vZ10 = [vhZ(hlid);vhZ(clid)];
vZ20 = [vhZ(hhid);vhZ(chid)];
vF10 = [vhF(hlid);vhF(clid)];
vF20 = [vhF(hhid);vhF(chid)];
vP10 = [vhP(hlid);vhP(clid)];
vP20 = [vhP(hhid);vhP(chid)];

mF10 = fitlm(vZ10,vF10,'linear');
mP10 = fitlm(vZ10,vP10,'linear');
mF20 = fitlm(vZ20,vF20,'linear');
mP20 = fitlm(vZ20,vP20,'linear');

zoo = [0:0.25:3]';
pF10 = predict(mF10,zoo);
pP10 = predict(mP10,zoo);
pF20 = predict(mF20,zoo);
pP20 = predict(mP20,zoo);


epred(:,1,j) = pF10;
epred(:,2,j) = pP10;
epred(:,3,j) = pF20;
epred(:,4,j) = pP20;

%%
% cm21=[1 0.5 0;...   %orange
%     0.5 0.5 0;... %tan/army
%     0 0.7 0;...   %g
%     0 1 1;...     %c
%     0 0 0.75;...  %b
%     0.5 0 1;...   %purple
%     1 0 1;...     %m
%     1 0 0;...     %r
%     0.5 0 0;...   %maroon
%     0.75 0.75 0.75;... %lt grey
%     0.5 0.5 0.5;...    %med grey
%     49/255 79/255 79/255;... %dk grey
%     0 0 0;...      %black
%     1 1 0;...      %yellow
%     127/255 255/255 0;... %lime green
%     0 0.5 0;...    %dk green
%     0/255 206/255 209/255;... %turq
%     0 0.5 0.75;...   %med blue
%     188/255 143/255 143/255;... %rosy brown
%     255/255 192/255 203/255;... %pink
%     255/255 160/255 122/255]; %peach

figure(2)
clf
subplot(2,2,1)
plot(vZ10,vP10,'o','color',[0/255 240/255 240/255]); hold on;
plot(vZ10,vF10,'x','color',[255/255 192/255 203/255]); hold on;
plot(zoo,pP10,'b','LineWidth',2); hold on;
plot(zoo,pF10,'r','LineWidth',2); hold on;
title('10^oC')
xlabel('zooplankton production (g m^-^2 d^-^1)')
ylabel('fish production (g m^-^2 d^-^1)')
axis([0 3 0 0.4])
legend('P','F')
legend('location','northwest')

subplot(2,2,2)
plot(vZ20,vP20,'o','color',[0/255 240/255 240/255]); hold on;
plot(vZ20,vF20,'x','color',[255/255 192/255 203/255]); hold on;
plot(zoo,pP20,'b','LineWidth',2); hold on;
plot(zoo,pF20,'r','LineWidth',2); hold on;
title('25^oC')
xlabel('zooplankton production (g m^-^2 d^-^1)')
ylabel('fish production (g m^-^2 d^-^1)')
axis([0 3 0 0.4])
%print('-dpng',[ppath 'Hist_Fore_',harv,'_FPprod_zProd_temp2_scatter_',simname,'.png'])

end
%%
save([nfile 'FPprod_zProd_temp2_hist_fore_fits_ens.mat'],'epred');


