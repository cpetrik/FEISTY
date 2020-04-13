% Calculate production / scope for growth
% (assim*Ingest - Metab)
% to see if emergent pattern of best and worst psets

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/',...
    'Matlab_New_sizes/param_ensemble/',...
    'Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/full_runs/'];

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath = ['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

%% Ensemble
dpath = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_enc-k063_met-k086_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
epath = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/Dc_enc-k063_cmax20-b250-k063_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/';

load([epath 'LHS_param6_mid6_kt3_bestAIC_params_Fupneg_mult10_Pneg2_mult3_reduced.mat'],...
    'red_params','ptext');
load([epath 'Hist_Fore_',harv,'_ensem_mid6_temp3_pset_VarMaxMinDiffSims_prod.mat']);

assim = repmat(pstats(:,2),1,3);
bpow  = repmat(pstats(:,3),1,3);
benc  = repmat(pstats(:,4),1,3);
amet  = repmat(pstats(:,5),1,3);
gam   = repmat(pstats(:,6),1,3);
kt    = repmat(pstats(:,7),1,3);

texp = {'xA','nA','xF','nF','xP','nP','xD','nD','xB','nB'};

%% Physiol
%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);
wgt=[M_s, M_m, M_l];
m = repmat(wgt,10,1);

Temps = [0,10,20,30];
temp  = 15;
Preys = 0:0.05:5;
prey = 2;
frac = 1;
pref = 1;

% Encounter rate %Me
%A = (exp(0.063*(temp-10.0)) .* 70 .* m.^(-0.20)) ./365.0;
A = (exp(0.063.*(temp-10.0)) .* gam .* m.^(-benc)) ./365.0;
enc = prey.*A.*frac.*pref;

% Max ingestion %Me
cmax = (exp(0.063.*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;
%cmax = (exp(kc*(temp-10.0)) .* h .* wgt^(-bcmx)) ./365.0;

% Ingestion
%ENC = sum(enc,2); % total biomass encountered
%con = cmax .* enc(:,1) ./ (cmax + ENC); % Type II
%type II fn
ENC = enc;
%con = (cmax .* enc) ./ (cmax + ENC);
con = cmax;
ing = assim .* con;

% Metabolism
met = (exp(kt*(temp-10.0)) .* amet .* m.^(-bpow)) ./365.0;

% scope for growth
sfg = ing - met;

%scale with cmax
flev = sfg ./ cmax;
smet = met ./ cmax;

%%
% figure
% subplot(2,2,1)
% bar(sfg(:,2))
% title('M')
% 
% subplot(2,2,2)
% bar(sfg(:,3))
% title('L')

figure(2)
subplot(2,2,1)
bar(flev(:,2))
ylim([0 0.6])
set(gca,'XTickLabel',texp)
ylabel('M nu scaled')
text(13,0.67,'con = cmax, T=15C','HorizontalAlignment','center')

subplot(2,2,2)
bar(flev(:,3))
ylim([0 0.6])
set(gca,'XTickLabel',texp)
ylabel('L nu scaled')

subplot(2,2,3)
bar(smet(:,2))
ylim([0 0.6])
set(gca,'XTickLabel',texp)
ylabel('M met scaled')

subplot(2,2,4)
bar(smet(:,3))
ylim([0 0.6])
set(gca,'XTickLabel',texp)
ylabel('L met scaled')
print('-dpng',[ppath 'Scaled_nu_met_MaxMinDiffSims_prod.png'])

%% Construct temp & prey combinations for M and L fishes only
Temps = 0:5:35;
%Preys = 0.2:0.2:5;

spath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' cfile '/'];
load([spath 'Prey_quant_Fore_Hist_All_fish03_' cfile '.mat'])
Preys = q(:,12);

% assim = repmat(pstats(:,2),1,length(Preys));
% bpow  = repmat(pstats(:,3),1,length(Preys));
% benc  = repmat(pstats(:,4),1,length(Preys));
% amet  = repmat(pstats(:,5),1,length(Preys));
% gam   = repmat(pstats(:,6),1,length(Preys));
% kt    = repmat(pstats(:,7),1,length(Preys));

assim = repmat(pstats(:,2),1,length(Temps));
bpow  = repmat(pstats(:,3),1,length(Temps));
benc  = repmat(pstats(:,4),1,length(Temps));
amet  = repmat(pstats(:,5),1,length(Temps));
gam   = repmat(pstats(:,6),1,length(Temps));
kt    = repmat(pstats(:,7),1,length(Temps));

% T = repmat(Temps,length(Preys),1);
% P = repmat(Preys,1,length(Temps));
T = repmat(Temps,10,1);
%P = repmat(Preys,10,1);

%% loop prey conc
nT = length(Temps);
nP = length(Preys);
Msfg = nan*ones(length(kt),length(Temps)+1,length(Preys)+1);
Lsfg = nan*ones(length(kt),length(Temps)+1,length(Preys)+1);
for i=1:length(Preys)
    P = Preys(i);
    
    % Encounter rate %Me
    MA = (exp(0.063.*(T-10.0)) .* gam .* M_m.^(-benc)) ./365.0;
    Menc = P.*MA.*frac.*pref;
    
    % Max ingestion %Me
    Mcmax = (exp(0.063.*(T-10.0)) .* 20 .* M_m.^(-0.25)) ./365.0;
    
    % Ingestion
    %type II fn
    MENC = Menc;
    Mcon = (Mcmax .* Menc) ./ (Mcmax + MENC);
    Ming = assim .* Mcon;
    
    % Metabolism
    Mmet = (exp(kt.*(T-10.0)) .* amet .* M_m.^(-bpow)) ./365.0;
    
    % scope for growth
    Msfg(:,1:nT,i) = (Ming - Mmet) ./ Mcmax;
    
    
    % Encounter rate
    LA = (exp(0.063.*(T-10.0)) .* gam .* M_l.^(-benc)) ./365.0;
    Lenc = P.*LA.*frac.*pref;
    
    % Max ingestion
    Lcmax = (exp(0.063.*(T-10.0)) .* 20 .* M_l.^(-0.25)) ./365.0;
    
    % Ingestion
    %type II fn
    LENC = Lenc;
    Lcon = (Lcmax .* Lenc) ./ (Lcmax + LENC);
    Ling = assim .* Lcon;
    
    % Metabolism
    Lmet = (exp(kt.*(T-10.0)) .* amet .* M_l.^(-bpow)) ./365.0;
    
    % scope for growth
    Lsfg(:,1:nT,i) = (Ling - Lmet) ./ Lcmax;
end

%% Pcolor plots
ptx = {'0','0.15','0.34','1.6','3.4','7.0'};

% All
figure(1)
subplot(2,2,1)
pcolor(squeeze(Msfg(1,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M max All')

subplot(2,2,2)
pcolor(squeeze(Msfg(2,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M min All')

subplot(2,2,3)
pcolor(squeeze(Lsfg(1,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L max All')

subplot(2,2,4)
pcolor(squeeze(Lsfg(2,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L min All')
print('-dpng',[ppath 'Scaled_nu_MaxMinDiffSims_prod_All.png'])

%% Forage
figure(2)
subplot(2,2,1)
pcolor(squeeze(Msfg(3,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M max F')

subplot(2,2,2)
pcolor(squeeze(Msfg(4,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M min F')

subplot(2,2,3)
pcolor(squeeze(Lsfg(3,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L max F')

subplot(2,2,4)
pcolor(squeeze(Lsfg(4,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L min F')
print('-dpng',[ppath 'Scaled_nu_MaxMinDiffSims_prod_F.png'])

%% Large pelagic
figure(3)
subplot(2,2,1)
pcolor(squeeze(Msfg(5,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M max P')

subplot(2,2,2)
pcolor(squeeze(Msfg(6,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M min P')

subplot(2,2,3)
pcolor(squeeze(Lsfg(5,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L max P')

subplot(2,2,4)
pcolor(squeeze(Lsfg(6,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L min P')
print('-dpng',[ppath 'Scaled_nu_MaxMinDiffSims_prod_P.png'])

%% Demersal
figure(4)
subplot(2,2,1)
pcolor(squeeze(Msfg(7,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M max D')

subplot(2,2,2)
pcolor(squeeze(Msfg(8,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M min D')

subplot(2,2,3)
pcolor(squeeze(Lsfg(7,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L max D')

subplot(2,2,4)
pcolor(squeeze(Lsfg(8,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L min D')
print('-dpng',[ppath 'Scaled_nu_MaxMinDiffSims_prod_D.png'])

%% Benthos
figure(5)
subplot(2,2,1)
pcolor(squeeze(Msfg(9,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M max B')

subplot(2,2,2)
pcolor(squeeze(Msfg(10,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('M min B')

subplot(2,2,3)
pcolor(squeeze(Lsfg(9,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L max B')

subplot(2,2,4)
pcolor(squeeze(Lsfg(10,:,:)))
shading flat
colorbar
cmocean('thermal')
xlabel('prey')
ylabel('temp')
set(gca,'YTick',1:8,'YTickLabel',Temps,'XTick',1:6,'XTickLabel',ptx)
caxis([0 0.6])
title('L min B')
print('-dpng',[ppath 'Scaled_nu_MaxMinDiffSims_prod_B.png'])

