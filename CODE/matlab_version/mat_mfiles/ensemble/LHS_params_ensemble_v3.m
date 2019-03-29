% Create 100 sets of parameters for FEISTY runs using
% Latin Hypercube Sampling
% reduced parameter ranges

clear all
close all

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/Dc_D075_Sm025_nmort1_noCC_RE00100/';

% ptext = {'frate','Lambda','K_a','amet','h','gam','kc','ke','kt','bpow','benc',...
%     'bcmx','bent_eff','J','A'};
% plow = [0.15;   0.56;   0.5;    2;  20;  20;  0.0302;0.0302;0.0302;0.1;   0.1;...
%     0.1;    0.05;  0.5;0.5];
% phi = [0.5;   0.84;     1;      8;  200;  200;  0.1208;0.1208;0.1208;0.3;0.3;...
%     0.3;    0.15;  1; 1];

ptext = {'Lambda','amet','h','gam','bpow','benc','bcmx','A'};
plow  = [0.56;       2;   20;  20;   0.1;   0.1;   0.1;  0.5];
phi   = [0.84;       8;  200; 200;   0.3;   0.3;   0.3;  1];

X=lhsdesign(200,length(ptext),'criterion','correlation');

%% orig set
plow = plow';
phi = phi';
pL = repmat(plow,length(X),1);
pH = repmat(phi,length(X),1);
fx = ((pH-pL) .* X) + pL;
%%
figure
for i=1:8
    subplot(3,3,i)
    hist(fx(:,i))
    title(ptext{i})
end
print('-dpng',[pp 'LHS_coverage_all.png'])

%% require b_met < b_cmax
id=find(X(:,5)<X(:,7));
Y = X(id,:);


%% set actual parameter values using ranges
pL = repmat(plow,length(id),1);
pH = repmat(phi,length(id),1);
fx = ((pH-pL) .* Y) + pL;

fmin = min(fx);
fmax = max(fx);

%%
figure
for i=1:8
    subplot(3,3,i)
    hist(fx(:,i))
    title(ptext{i})
end
print('-dpng',[pp 'LHS_coverage_bM_bC_force.png'])

%% Save
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
save([nfile 'LHS_param15_100vals_v3.mat'],'fx','ptext','plow','phi');









