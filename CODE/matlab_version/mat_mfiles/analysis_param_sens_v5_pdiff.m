% POEM output at all locations
% Only last 12 months of 150 years saved 
% Biomass in temperate & polar vs. subtropics and tropics
% instead of mean latitude

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];

load([nfile 'Climatol_All_fish03_means_param_sens_v3.mat'])
%ptext{29}='kap75';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%% Outputs
F = SF + MF;
P = SP + MP + LP;
D = SD + MD + LD;
All = F + P + D;
% Percent change in biomass of each group
v1 = (nanmean(F) - nanmean(F(:,1))) ./ ((nanmean(F) + nanmean(F(:,1)))/2);
v2 = (nanmean(P) - nanmean(P(:,1))) ./ ((nanmean(P) + nanmean(P(:,1)))/2);
v3 = (nanmean(D) - nanmean(D(:,1))) ./ ((nanmean(D) + nanmean(D(:,1)))/2);
% Total biomass
v4 = (nanmean(All) - nanmean(All(:,1))) ./ ...
    ((nanmean(All) + nanmean(All(:,1)))/2);

%% Biomass in tropics vs. out
vlat = lat(ID);
trop = (abs(vlat)<=30);
eq = nanmean(All(trop==1,:));
pol = nanmean(All(trop==0,:));
v5 = (eq - eq(:,1)) ./ ((eq + eq(:,1))/2);
v6 = (pol - pol(:,1)) ./ ((pol + pol(:,1))/2);

% Combined
v = sqrt(v1.^2 + v2.^2 + v3.^2 + v5.^2 + v6.^2);

%% Bar graph
figure(1)
subplot(3,2,1)
bar((v1))
title('F')
subplot(3,2,3)
bar((v2))
title('P')
subplot(3,2,5)
bar((v3))
title('D')
ylabel('percent difference')
subplot(3,2,2)
bar((v5))
title('Tropics')
subplot(3,2,4)
bar((v6))
title('Temp-Polar')
subplot(3,2,6)
bar((v))
title('total magnitude')
%%
figure(2)
subplot(3,2,1)
bar(v1(3:end))
subplot(3,2,2)
bar(v2(3:end))
subplot(3,2,3)
bar(v3(3:end))
subplot(3,2,4)
bar(v5(3:end))
subplot(3,2,5)
bar(v6(3:end))
subplot(3,2,6)
bar(v(3:end))

%% Orthogonality
vall(1,:) = v1;
vall(2,:) = v2;
vall(3,:) = v3;
vall(4,:) = v5;
vall(5,:) = v6;
vall(6,:) = v4;

ortho = NaN*ones(length(v));
for n=1:length(v)
    for m=1:length(v)
        ortho(m,n) = 1 - ((vall(:,m).' * vall(:,n)) / (v(m) * v(n)));
    end
end

%% Cluster analysis (in R)
vtab = table(v1',v2',v3',v6',v5',v4','VariableNames',{'F','P','D','Trop',...
    'Temp','All'},'RowNames',ptext');
writetable(vtab,[nfile 'Climatol_All_fish03_param_sens_v5_vecs_pdiff.csv'],...
    'WriteRowNames',true)




