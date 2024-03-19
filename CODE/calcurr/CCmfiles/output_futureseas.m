% FEISTY output for Ensemble analysis

% 1. Total biomass of sardine
% 2. Growth/maturation from one size class to the next = rec
% 3. Energy available for growth scaled by the biomass = prod
% 4. Population prod = biomass growing/maturing from the juvenile size 
% class to the adult size class at time=T divided by the biomass of adults 
% at time=(T-1)
% 5. Catch N of 42nd parallel 
% 6. Catch S of 42nd parallel 
% 7. Annual predation mortality rate
% 8. Ecosystem prod = Forage biomass / Zooplankton biomass

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_CC80_RE00100';
harv = 'All_fishobs';
esm1 = {'IPSL','GFDL','HAD'};
esm2 = {'ipsl','gfdl','hadley'};

vers = 'IPSL';
fpath=['/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/',vers,'/'];

%%
load([fpath 'Means_Hist_' vers '_' harv '_' cfile '.mat'])
load([fpath 'Means_Hist_' vers '_' harv '_prod_' cfile '.mat'])
load([fpath 'Means_Hist_' vers '_' harv '_rec_' cfile '.mat'])
load([fpath 'Means_Hist_' vers '_' harv '_mort_' cfile '.mat'])

%% Need to do area-weighted means & sums

%% means or sums
nt = length(time);
% mean prod
f_tmprod = (sf_tmprod + mf_tmprod)/2;
% pop prod
f_pp = NaN*ones(nt,1);
f_pp(2:end) = mf_tmrec(2:end) ./ mf_tmean(1:(end-1));
% pred mort
f_tmpred = (sf_tmpred + mf_tmpred)/2;

%% Each year
a = 1:12:nt; % start of each yr
b = 12:12:nt; % end of each yr

aFrec = NaN*ones((nt/12),1);
aFprod = aFrec;
aFpopprod = aFrec;
aFpred = aFrec;

for i = 1:(nt/12)
    
    aFrec(i) = mean(mf_tmrec(a(i):b(i)),2,"omitnan");
    aFprod(i) = mean(f_tmprod(a(i):b(i)),2,"omitnan");
    aFpred(i) = mean(mf_tmpred(a(i):b(i)),2,"omitnan");
    aFpopprod(i) = mean(f_pp(a(i):b(i)),1,"omitnan");
    
end

%% Create csv table
modyrs = 1980:2010;





