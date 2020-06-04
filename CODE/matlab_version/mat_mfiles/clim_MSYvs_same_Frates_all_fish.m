% FEISTY F vs. MMSSY
% Same F rate applied to all fish
% Climatology

clear all
close all

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];

ppath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/';

load([dp 'Climatol_same_Frates_all_fish.mat']);

%% total catches
lme_Amcatch = lme_Fmcatch + lme_Pmcatch + lme_Dmcatch;

Acatch = sum(MFc+MPc+MDc+LPc+LDc);
Fcatch = sum(MFc);
Pcatch = sum(MPc+LPc);
Dcatch = sum(MDc+LDc);

%% Plots of catch vs. F
% Global
figure(1)
subplot(2,2,1)
plot(Fish,Fcatch,'.r','MarkerSize',10)
ylabel('catch (g m^-^2 y^-^1)')
title('Forage')

subplot(2,2,2)
plot(Fish,Pcatch,'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
plot(Fish,Dcatch,'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('F')
ylabel('catch (g m^-^2 y^-^1)')
title('Demersal')

subplot(2,2,4)
plot(Fish,Acatch,'.k','MarkerSize',10)
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_global_yield_sameF_all_fish.png'])

%% LME
nrows=11;
ncols=6;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);

%% all 
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Amcatch(i,:))*1e-12,'.k','MarkerSize',10); hold on;
        xlim([0 5])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_All.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Fmcatch(i,:))*1e-12,'.r','MarkerSize',10); hold on;
        xlim([0 5])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Pmcatch(i,:))*1e-12,'.b','MarkerSize',10); hold on;
        xlim([0 5])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_P.png'])

%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Dmcatch(i,:))*1e-12,'.','color',[0 0.6 0],'MarkerSize',10); hold on;
        xlim([0 5])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_sameF_all_fish_D.png'])

%% Find MSY
fid = max(lme_Fmcatch,[],2);
pid = max(lme_Pmcatch,[],2);
did = max(lme_Dmcatch,[],2);
aid = max(lme_Amcatch,[],2);

Fmsy = NaN*ones(66,4);
for n=1:66
    Fmsy(n,1) = Fish(find(lme_Fmcatch(n,:)==fid(n)));
    Fmsy(n,2) = Fish(find(lme_Pmcatch(n,:)==pid(n)));
    Fmsy(n,3) = Fish(find(lme_Dmcatch(n,:)==did(n)));
    Fmsy(n,4) = Fish(find(lme_Amcatch(n,:)==aid(n)));
end

Fmsy5 = (Fmsy==5);

%% Plots of comp with SAUP
%Corr
figure(6)
subplot(2,2,1)
semilogx(Fish,r_all(2,:),'.r','MarkerSize',10)
ylabel('r')
title('Forage')

subplot(2,2,2)
semilogx(Fish,r_all(3,:),'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
semilogx(Fish,r_all(4,:),'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('F')
ylabel('r')
title('Demersal')

subplot(2,2,4)
semilogx(Fish,r_all(1,:),'.k','MarkerSize',10)
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_yield_sameF_all_fish_SAUPcomp_rcorr.png'])

%% RMSE
figure(7)
subplot(2,2,1)
semilogx(Fish,rmse_all(2,:),'.r','MarkerSize',10)
xlim([Fish(2) 5])
ylabel('RMSE')
title('Forage')

subplot(2,2,2)
semilogx(Fish,rmse_all(3,:),'.b','MarkerSize',10)
xlim([Fish(2) 5])
title('Large pelagic')

subplot(2,2,3)
semilogx(Fish,rmse_all(4,:),'.','color',[0 0.6 0],'MarkerSize',10)
xlim([Fish(2) 5])
xlabel('F')
ylabel('RMSE')
title('Demersal')

subplot(2,2,4)
%plot(Fish,rmse_all(1,:),'.k','MarkerSize',10)
semilogx(Fish,rmse_all(1,:),'.k','MarkerSize',10)
xlim([Fish(2) 5])
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_yield_sameF_all_fish_SAUPcomp_rmse.png'])

%% SS
figure(8)
subplot(2,2,1)
semilogx(Fish,ss_all(2,:),'.r','MarkerSize',10)
ylabel('Sum of squares')
title('Forage')

subplot(2,2,2)
semilogx(Fish,ss_all(3,:),'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
semilogx(Fish,ss_all(4,:),'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('F')
ylabel('Sum of squares')
title('Demersal')

subplot(2,2,4)
semilogx(Fish,ss_all(1,:),'.k','MarkerSize',10)
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_yield_sameF_all_fish_SAUPcomp_ss.png'])

