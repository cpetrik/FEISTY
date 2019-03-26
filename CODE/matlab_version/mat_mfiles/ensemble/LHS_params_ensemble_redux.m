% Create 100 sets of parameters for FEISTY runs using
% Latin Hypercube Sampling
% reduced parameter ranges

clear all
close all

ptext = {'frate','Lambda','K_a','amet','h','gam','kc','ke','kt','bpow','benc',...
    'bcmx','bent_eff','J','A'};
plow = [0.15;   0.56;   0.5;    2;  20;  20;  0.0302;0.0302;0.0302;0.1;   0.1;...
    0.1;    0.05;  0.5;0.5];
phi = [0.5;   0.84;     1;      8;  200;  200;  0.1208;0.1208;0.1208;0.3;0.3;...
    0.3;    0.15;  1; 1];

X=lhsdesign(100,length(ptext),'criterion','correlation');

%% check param coverage
r=corrcoef(X);
figure
pcolor(r)
caxis([-0.2 0.2])
colorbar

figure
for i=1:15
    subplot(4,4,i)
    hist(X(:,i))
end

%% set actual parameter values using ranges
plow = plow';
phi = phi';
pL = repmat(plow,100,1);
pH = repmat(phi,100,1);
fx = ((pH-pL) .* X) + pL;

fmin = min(fx);
fmax = max(fx);

%% Save
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
save([nfile 'LHS_param15_100vals_redux.mat'],'fx','ptext','plow','phi');









