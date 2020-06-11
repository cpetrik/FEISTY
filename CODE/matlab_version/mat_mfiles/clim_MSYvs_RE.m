% FEISTY F vs. MMSSY
% Same F rate applied to all fish
% Climatology

clear all
close all

dp = '/Volumes/FEISTY/NC/Matlab_new_size/';

ppath='/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/bio_rates/';

load([dp 'Climatol_All_fish03_same_F_RE_sens.mat']);

%% total catches
lme_Amcatch = lme_Fmcatch + lme_Pmcatch + lme_Dmcatch;

Acatch = sum(MFc+MPc+MDc+LPc+LDc);
Fcatch = sum(MFc);
Pcatch = sum(MPc+LPc);
Dcatch = sum(MDc+LDc);

%% Plots of catch vs. RE
% Global
figure(1)
subplot(2,2,1)
semilogx(reff,Fcatch,'.r','MarkerSize',10)
ylabel('catch (g m^-^2 y^-^1)')
title('Forage')

subplot(2,2,2)
semilogx(reff,Pcatch,'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
semilogx(reff,Dcatch,'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('F')
ylabel('catch (g m^-^2 y^-^1)')
title('Demersal')

subplot(2,2,4)
semilogx(reff,Acatch,'.k','MarkerSize',10)
xlabel('RE')
title('All')
print('-dpng',[ppath 'Climatol_All_fish03_global_yield_vsRE.png'])

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
        semilogx(reff,(lme_Amcatch(i,:))*1e-12,'.k','MarkerSize',10); hold on;
        %xlim([0 5])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_All_fish03_LME_yield_vsRE_All.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        semilogx(reff,(lme_Fmcatch(i,:))*1e-12,'.r','MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_All_fish03_LME_yield_vsRE_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        semilogx(reff,(lme_Pmcatch(i,:))*1e-12,'.b','MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_All_fish03_LME_yield_vsRE_P.png'])

%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        semilogx(reff,(lme_Dmcatch(i,:))*1e-12,'.','color',[0 0.6 0],'MarkerSize',10); hold on;
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_All_fish03_LME_yield_vsRE_D.png'])

%% Find MSY
fid = max(lme_Fmcatch,[],2);
pid = max(lme_Pmcatch,[],2);
did = max(lme_Dmcatch,[],2);
aid = max(lme_Amcatch,[],2);

Fmsy = NaN*ones(66,4);
for n=1:66
    Fmsy(n,1) = reff(find(lme_Fmcatch(n,:)==fid(n)));
    Fmsy(n,2) = reff(find(lme_Pmcatch(n,:)==pid(n)));
    Fmsy(n,3) = reff(find(lme_Dmcatch(n,:)==did(n)));
    Fmsy(n,4) = reff(find(lme_Amcatch(n,:)==aid(n)));
end

Fmsy5 = (Fmsy==reff(end));
sum(Fmsy5)
mean(Fmsy)

%% Plots of comp with SAUP
%Corr
figure(6)
subplot(2,2,1)
semilogx(reff,r_all(2,:),'.r','MarkerSize',10)
ylabel('r')
title('Forage')

subplot(2,2,2)
semilogx(reff,r_all(3,:),'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
semilogx(reff,r_all(4,:),'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('RE')
ylabel('r')
title('Demersal')

subplot(2,2,4)
semilogx(reff,r_all(1,:),'.k','MarkerSize',10)
xlabel('RE')
title('All')
print('-dpng',[ppath 'Climatol_All_fish03_yield_vsRE_SAUPcomp_rcorr.png'])

%% RMSE
figure(7)
subplot(2,2,1)
semilogx(reff,rmse_all(2,:),'.r','MarkerSize',10)
ylabel('RMSE')
title('Forage')

subplot(2,2,2)
semilogx(reff,rmse_all(3,:),'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
semilogx(reff,rmse_all(4,:),'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('RE')
ylabel('RMSE')
title('Demersal')

subplot(2,2,4)
%plot(reff,rmse_all(1,:),'.k','MarkerSize',10)
semilogx(reff,rmse_all(1,:),'.k','MarkerSize',10)
xlabel('RE')
title('All')
print('-dpng',[ppath 'Climatol_All_fish03_yield_vsRE_SAUPcomp_rmse.png'])

%% SS
figure(8)
subplot(2,2,1)
semilogx(reff,ss_all(2,:),'.r','MarkerSize',10)
ylabel('Sum of squares')
title('Forage')

subplot(2,2,2)
semilogx(reff,ss_all(3,:),'.b','MarkerSize',10)
title('Large pelagic')

subplot(2,2,3)
semilogx(reff,ss_all(4,:),'.','color',[0 0.6 0],'MarkerSize',10)
xlabel('RE')
ylabel('Sum of squares')
title('Demersal')

subplot(2,2,4)
semilogx(reff,ss_all(1,:),'.k','MarkerSize',10)
xlabel('RE')
title('All')
print('-dpng',[ppath 'Climatol_All_fish03_yield_vsRE_SAUPcomp_ss.png'])

