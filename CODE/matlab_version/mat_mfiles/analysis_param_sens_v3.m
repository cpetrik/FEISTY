% POEM output at all locations
% Only last 12 months of 150 years saved 

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
% Relative biomass of each group
v1 = nanmean(log10(F) - log10(F(:,1)));
v2 = nanmean(log10(P) - log10(P(:,1)));
v3 = nanmean(log10(D) - log10(D(:,1)));
% Total biomass
v4 = nanmean(log10(All) - log10(All(:,1)));
% Weighted latitudinal center of biomass
vlat = lat(ID);
nhem = find(vlat>0);
shem = find(vlat<0);
totb = repmat(nansum(All),length(ID),1);
nh = All(nhem,:) .* vlat(nhem) ./ totb(nhem,:);
NH = nansum(nh);
sh = All(shem,:) .* vlat(shem) ./ totb(shem,:);
SH = nansum(sh);
latc = (NH-SH)/2;
v5 = log10(latc) - log10(latc(:,1));

% Combined
v = sqrt(v1.^2 + v2.^2 + v3.^2 + v4.^2 + v5.^2);

%% Bar graph
figure(1)
subplot(3,2,1)
bar(v1)
subplot(3,2,2)
bar(v2)
subplot(3,2,3)
bar(v3)
subplot(3,2,4)
bar(v4)
subplot(3,2,5)
bar(v5)
subplot(3,2,6)
bar(v)

figure(2)
subplot(3,2,1)
bar(v1(3:end))
subplot(3,2,2)
bar(v2(3:end))
subplot(3,2,3)
bar(v3(3:end))
subplot(3,2,4)
bar(v4(3:end))
subplot(3,2,5)
bar(v5(3:end))
subplot(3,2,6)
bar(v(3:end))

%% Orthogonality
vall(1,:) = v1;
vall(2,:) = v2;
vall(3,:) = v3;
vall(4,:) = v4;
vall(5,:) = v5;

ortho = NaN*ones(length(v));
for n=1:length(v)
    for m=1:length(v)
        ortho(m,n) = 1 - ((vall(:,m).' * vall(:,n)) / (v(m) * v(n)));
    end
end

%% Cluster analysis (in R)
vtab = table(v1',v2',v3',v4',v5','VariableNames',{'F','P','D','All','lat'},...
    'RowNames',ptext');
writetable(vtab,[nfile 'Climatol_All_fish03_param_sens_v3_vecs.csv'],...
    'WriteRowNames',true)




