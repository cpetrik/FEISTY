% Create full sets of parameters for FEISTY runs using
% permn to get all 2^5 combos of hi/low for 5 parameters

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

ptext = {'Lambda','bpow','benc','amet','gam'};
plow  = [0.6;     0.15;   0.15;  3;   50];
phi   = [0.75;    0.25;   0.25;  6;   100];


%% set actual parameter values using ranges
X=permn([0 1],5);
pL = repmat(plow,length(X),1);
pH = repmat(phi,length(X),1);
fx = ((pH-pL) .* X) + pL;

%%
figure
for i=1:5
    subplot(2,3,i)
    hist(fx(:,i))
    title(ptext{i})
end
%print('-dpng',[pp 'LHS_coverage_all.png'])

%% Save
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
save([nfile 'LHS_param5_hi_low.mat'],'fx','ptext','plow','phi');


