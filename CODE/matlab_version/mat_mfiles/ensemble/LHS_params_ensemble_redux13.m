% Create 100 sets of parameters for FEISTY runs using
% Latin Hypercube Sampling
% reduced parameter ranges using results of 5 most sens params hi/low

clear all
close all

ptext = {'Lambda','K_a','amet','h','gam','kc','ke','kt','bpow','benc',...
    'bcmx','bent_eff','A'};
plow = [0.6;   0.5;    2;  10;  35;  0.030;0.030;0.030; 0.1;    0.15;...
    0.2;    0.05;  0.5];
phi = [0.8;     1;     6;  40;  140; 0.121;0.121;0.121; 0.199;   0.25;...
    0.3;    0.15;  1];

X=lhsdesign(100,length(ptext),'criterion','correlation');

%% check param coverage
r=corrcoef(X);
figure
pcolor(r)
caxis([-0.2 0.2])
colorbar

figure
for i=1:13
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
%%
figure
for i=1:13
    subplot(4,4,i)
    hist(fx(:,i))
    title(ptext{i})
end

%% Save
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
save([nfile 'LHS_param13_100vals.mat'],'fx','ptext','plow','phi');









