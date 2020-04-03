% Create full sets of parameters for FEISTY runs using
% permn to get all 3^6 combos of hi/mid/low for 6 parameters
% kts of 0.063 and 0.108 did not have good AICs (b/c of coexist)
% try k = 0.063, mid, and 0.0955 for all rates (ignoring coexist)

clear all
close all

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/'];

ptext2 = {'Lambda','bpow','benc','amet','gam','kt'};
plow2  = [0.6,     0.15,   0.15,  3,   50,  0.063];
phi2   = [0.75,    0.20,   0.25,  5,   100, 0.0955];


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
% plow3  = [0.6,     0.15,   0.15,  3,   50];
% phi3   = [0.75,    0.20,   0.25,  5,   100];
% X3=permn([0 0.5 1],5);
% 
% % set actual parameter values using ranges
% pL3 = repmat(plow3,length(X3),1);
% pH3 = repmat(phi3,length(X3),1);
% fz3 = ((pH3-pL3) .* X3) + pL3;
% fz3(:,6) = 0.0955;
% 
% fx_all = [fz; fz3];
% id = find(fx_all(:,6)>0.09);
% fx = fx_all(id,:);

fx = fz;

%% Save
save([nfile 'LHS_param6_mid6_samek.mat']);


