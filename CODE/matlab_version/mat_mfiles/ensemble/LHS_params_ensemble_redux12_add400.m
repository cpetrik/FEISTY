% Add 400 sets to LHS 100 sets of parameters for FEISTY runs 
% reduced parameter ranges using results of 5 most sens params hi/low
% reduced again because first round was not enough

clear all
close all

ptext = {'Lambda','K_a','amet','h','gam','kce','kt','bpow','benc',...
    'bcmx','bent_eff','A'};
plow = [0.6;   0.5;    2;  10;  50;  0.030;0.030; 0.1;    0.15;...
    0.2;    0.05;  0.5];
phi = [0.75;  0.95;    6;  40; 150; 0.121;0.121; 0.199;   0.25;...
    0.3;    0.15;  0.8];

%X=lhsdesign(100,length(ptext),'criterion','correlation');
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/param_ensemble/'];
load([nfile 'LHS_param12_100vals.mat'],'fx');

plow = plow';
phi = phi';
pL = repmat(plow,100,1);
pH = repmat(phi,100,1);
%fx = ((pH-pL) .* X) + pL;
%fx - pL = (pH-pL) .* X
X = (fx - pL) ./ (pH-pL);

%%
figure
for i=1:12
    subplot(4,3,i)
    hist(X(:,i))
end

%%
xF=lhsaugment(X,400);

%%
figure
for i=1:12
    subplot(4,3,i)
    hist(xF(:,i))
end

%% set actual parameter values using ranges
pL = repmat(plow,500,1);
pH = repmat(phi,500,1);
fx2 = ((pH-pL) .* xF) + pL;

fmin = min(fx2);
fmax = max(fx2);
%%
figure
for i=1:12
    subplot(4,3,i)
    hist(fx2(:,i))
    title(ptext{i})
end

%% Save
save([nfile 'LHS_param12_500vals.mat'],'fx2','ptext','plow','phi');
