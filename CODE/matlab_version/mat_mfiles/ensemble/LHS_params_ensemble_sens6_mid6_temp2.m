% Create full sets of parameters for FEISTY runs using
% permn to get all 3^6 combos of hi/mid/low for 6 parameters
% kts of 0.063 and 0.108 did not have good AICs

clear all
close all

%% old one
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
load([nfile 'LHS_param5_mid5.mat']);

%%
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

ptext2 = {'Lambda','bpow','benc','amet','gam','kt'};
plow2  = [0.6,     0.15,   0.15,  3,   50, 0.0730];
phi2   = [0.75,    0.20,   0.25,  5,   100, 0.0980];


%% set actual parameter values using ranges
X=permn([0 0.5 1],6);
imagesc(X)

pL = repmat(plow2,length(X),1);
pH = repmat(phi2,length(X),1);
fz = ((pH-pL) .* X) + pL;

%%
figure
for i=1:6
    subplot(2,3,i)
    hist(fz(:,i))
    title(ptext2{i})
end
%print('-dpng',[pp 'LHS_coverage_all.png'])

%%
id1 = find(fz(:,6)<0.08);
id2 = find(fz(:,6)>0.09);
id = [id1;id2];
fx_all = fz;
fx = fz(id,:);

%% Save
save([nfile 'LHS_param6_mid6_temp2.mat']);


