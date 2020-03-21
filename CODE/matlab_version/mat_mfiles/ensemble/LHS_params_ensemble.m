% Create 100 sets of parameters for FEISTY runs using
% Latin Hypercube Sampling

clear all
close all

nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/'];


ptext = {'frate','Lambda','K_a','amet','h','gam','kc','ke','kt','bpow','benc',...
    'bcmx','bent_eff','J','A'};
plow = [0.15;   0.56;   0.25;   2;  5;  5;  0.0302;0.0302;0.0302;0.1;   0.1;...
    0.1;    0.0375;  0.5;0.5];
phi = [0.6;   0.84;     1;      8;  500;  500;  0.1208;0.1208;0.1208;0.32;0.32;...
    0.32;    0.15;  1; 1];

X=lhsdesign(100,length(ptext),'criterion','correlation');

ptab(:,1) = plow;
ptab(:,2) = phi;
pT = array2table(ptab,'VariableNames',{'Low','High'},'RowNames',ptext);
writetable(pT,[nfile 'LHS_param15_high_low.csv'],'WriteRowNames',true)


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
save([nfile 'LHS_param15_100vals.mat'],'fx','ptext','plow','phi');









