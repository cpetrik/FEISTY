% Create full sets of F combinations for 3 types
% Instead of 3 loops

clear all
close all

%%
nfile = ['/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/'];

ptext2 = {'FF','PF','DF'};

% set actual parameter values using ranges
X=permn([0:0.1:2],3);
imagesc(X)


%% Save
save([nfile 'LHS_diffF.mat']);


