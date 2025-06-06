% Create 100 sets of parameters for FEISTY runs using
% Latin Hypercube Sampling
% reduced parameter ranges using results of 5 most sens params hi/low
% reduced again because first round was not enough

clear all
close all

ptext = {'Lambda','amet','gam','bpow','benc','bcmx','bent_eff','A'};
plow =    [0.6;    2;     50;  0.1;    0.15;  0.2;    0.05;  0.5];
phi =     [0.725;  5;    150;  0.199;  0.25;  0.3;    0.15;  0.75];

X=lhsdesign(100,length(ptext),'criterion','correlation');

%% check param coverage
r=corrcoef(X);
figure
pcolor(r)
caxis([-0.2 0.2])
colorbar

figure
for i=1:8
    subplot(3,3,i)
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
for i=1:8
    subplot(3,3,i)
    hist(fx(:,i))
    title(ptext{i})
end

%% Save
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
save([nfile 'LHS_param8_100vals.mat'],'fx','ptext','plow','phi');









