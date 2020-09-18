% FEISTY F vs. MMSSY
% Diff F rate applied to one types
% Other types held constant at 0.2 or 0.3
% Climatology

clear all
close all

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/Climatology/'];

ppath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/',...
    'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

load([dp 'Climatol_diff_Frates_P2.mat']);

%% total catches
lme_Amcatch = lme_Fmcatch + lme_Pmcatch + lme_Dmcatch;

Acatch = sum(MFc+MPc+MDc+LPc+LDc);
Fcatch = sum(MFc);
Pcatch = sum(MPc+LPc);
Dcatch = sum(MDc+LDc);

Fish = 0:0.1:2;

%% Plots of catch vs. F
% Global
figure(1)
subplot(2,2,1)
plot(Fish,Fcatch,'.r','MarkerSize',20)
ylabel('catch (g m^-^2 y^-^1)')
title('Forage = 0.2')

subplot(2,2,2)
plot(Fish,Pcatch,'.b','MarkerSize',20)
title('Large pelagic varied')

subplot(2,2,3)
plot(Fish,Dcatch,'.','color',[0 0.6 0],'MarkerSize',20)
xlabel('F')
ylabel('catch (g m^-^2 y^-^1)')
title('Demersal = 0.2')

subplot(2,2,4)
plot(Fish,Acatch,'.k','MarkerSize',20)
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_global_yield_diff_Frates_P2.png'])

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
        xlim([0 2])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_diff_Frates_P2_All.png'])

%% F
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Fmcatch(i,:))*1e-12,'.r','MarkerSize',10); hold on;
        xlim([0 2])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_diff_Frates_P2_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Pmcatch(i,:))*1e-12,'.b','MarkerSize',10); hold on;
        xlim([0 2])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_diff_Frates_P2_P.png'])

%% D
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Dmcatch(i,:))*1e-12,'.','color',[0 0.6 0],'MarkerSize',10); hold on;
        xlim([0 2])
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'Climatol_LME_yield_diff_Frates_P2_D.png'])

%% Find MSY
fid = max(lme_Fmcatch,[],2);
pid = max(lme_Pmcatch,[],2);
did = max(lme_Dmcatch,[],2);
aid = max(lme_Amcatch,[],2);

Fmsy = NaN*ones(66,4);
for n=1:66
    fid2 = Fish(find(lme_Fmcatch(n,:)==fid(n)));
    pid2 = Fish(find(lme_Pmcatch(n,:)==pid(n)));
    did2 = Fish(find(lme_Dmcatch(n,:)==did(n)));
    aid2 = Fish(find(lme_Amcatch(n,:)==aid(n)));
    Fmsy(n,1) = fid2(1);
    Fmsy(n,2) = pid2(1);
    Fmsy(n,3) = did2(1);
    Fmsy(n,4) = aid2(1);
end

Fmsy5 = (Fmsy==2);
sum(Fmsy5)

%% Plots of comp with SAUP
%Corr
figure(6)
subplot(2,2,1)
plot(Fish,r_all(2,:),'.r','MarkerSize',20)
ylabel('r')
title('Forage = 0.2')

subplot(2,2,2)
plot(Fish,r_all(3,:),'.b','MarkerSize',20)
title('Large pelagic varied')

subplot(2,2,3)
plot(Fish,r_all(4,:),'.','color',[0 0.6 0],'MarkerSize',20)
xlabel('F')
ylabel('r')
title('Demersal = 0.2')

subplot(2,2,4)
plot(Fish,r_all(1,:),'.k','MarkerSize',20)
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_yield_diff_Frates_P2_SAUPcomp_rcorr.png'])

%% RMSE
figure(7)
subplot(2,2,1)
plot(Fish,rmse_all(2,:),'.r','MarkerSize',20)
%xlim([Fish(2) 5])
xlim([0 2])
ylabel('RMSE')
title('Forage = 0.2')

subplot(2,2,2)
plot(Fish,rmse_all(3,:),'.b','MarkerSize',20)
xlim([0.1 2])
title('Large pelagic varied')

subplot(2,2,3)
plot(Fish,rmse_all(4,:),'.','color',[0 0.6 0],'MarkerSize',20)
xlim([0 2])
xlabel('F')
ylabel('RMSE')
title('Demersal = 0.2')

subplot(2,2,4)
%plot(Fish,rmse_all(1,:),'.k','MarkerSize',20)
plot(Fish,rmse_all(1,:),'.k','MarkerSize',20)
xlim([0.1 2])
xlabel('F')
title('All')
print('-dpng',[ppath 'Climatol_yield_diff_Frates_P2_SAUPcomp_rmse.png'])

%% Subset of LMEs
lmes = [1:9, 54:55, 65, 21:22, 59:60];
nrows=4;
ncols=4;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);

%% F
figure('Units','inches','Position',[1 3 6.5 6.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        L = lmes(i);
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Fmcatch(L,:))*1e-12,'.-r','MarkerSize',10); hold on;
        xlim([0 2])
        title(num2str(L))
        
    end
end
print('-dpng',[ppath 'Climatol_LME_select_yield_diff_Frates_P2_F.png'])

%% P
figure('Units','inches','Position',[1 3 6.5 6.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        L = lmes(i);
        subplot('position',pos(m,:,n))
        plot(Fish,(lme_Dmcatch(L,:))*1e-12,'.-','color',[0 0.6 0],'MarkerSize',10); hold on;
        xlim([0 2])
        title(num2str(L))
        set(gca,'XTickLabel',[])
        
    end
end
print('-dpng',[ppath 'Climatol_LME_select_yield_diff_Frates_P2_D.png'])

%% Individual plots, no axis labels for ppt
% F
figure(15);
for i = 1:length(lmes)
        L = lmes(i);
        plot(Fish,(lme_Fmcatch(L,:))*1e-12,'.-r',...
            'MarkerSize',20,'LineWidth',2); 
        xlim([0 2])
        set(gca,'YTickLabel',[])
        title(num2str(L))
        ylabel('F Catch')
        xlabel('Fishing rate on P (yr^-^1)')
        print('-dpng',[ppath 'LMEs/Climatol_LME',num2str(L),'_yield_diff_Frates_P2_F.png'])
end

%% D
figure(16);
for i = 1:length(lmes)
        L = lmes(i);
        plot(Fish,(lme_Dmcatch(L,:))*1e-12,'.-','color',[0 0.6 0],...
            'MarkerSize',20,'LineWidth',2); 
        xlim([0 2])
        set(gca,'YTickLabel',[])
        title(num2str(L))
        ylabel('D Catch')
        xlabel('Fishing rate on P (yr^-^1)')
        print('-dpng',[ppath 'LMEs/Climatol_LME',num2str(L),'_yield_diff_Frates_P2_D.png'])
end

