% Create full sets of parameters for FEISTY runs using
% permn to get all 2^5 combos of hi/low for 5 parameters
% add midpoints for Lambda, benc, gam

clear all
close all

nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

%       'Lambda','benc','gam'
plow  = [0.6,     0.15,  50];
phi   = [0.75,    0.25,  100];

%% set actual parameter values using ranges
X=permn([0 0.5 1],3);
imagesc(X)

pL = repmat(plow,length(X),1);
pH = repmat(phi,length(X),1);

%%
ptext = {'Lambda','bpow','benc','amet','gam'};
fx = NaN(27,5);
fx(:,2) = 0.175*ones(27,1);
fx(:,4) = 4*ones(27,1);

x = ((pH-pL) .* X) + pL;
fx(:,1) = x(:,1);
fx(:,3) = x(:,2);
fx(:,5) = x(:,3);

%%
figure
for i=1:5
    subplot(2,3,i)
    hist(fx(:,i))
    title(ptext{i})
end
%print('-dpng',[pp 'LHS_coverage_all.png'])

%% Save
save([nfile 'LHS_param5_hi_low_mid_v2.mat'],'fx','ptext','plow','phi');


