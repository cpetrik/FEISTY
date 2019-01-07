% Visualize output of POEM
% Spinup at one location
% 100 years, but only last year saved

clear all
close all

%data path
cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC050_RE1000';
dpath = ['/Volumes/GFDL/CSV/Matlab_big_size/' cfile '/'];
%figures path
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/';
fpath = [pp cfile '/'];

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);
R=cell(length(spots),3);
%%
for s=1:length(spots)
    %%
    close all
    loc = spots{s};
    lname = [loc '_'];
    
    %     SF = Spinup_Sml_f;
    %     SP = Spinup_Sml_p;
    %     SD = Spinup_Sml_d;
    %     MF = Spinup_Med_f;
    %     MP = Spinup_Med_p;
    %     MD = Spinup_Med_d;
    %     LP = Spinup_Lrg_p;
    %     LD = Spinup_Lrg_d;
    %     C  = Spinup_Cobalt;
    
    
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    C = csvread([dpath sname lname 'Cobalt.csv']);
    z(:,1) = C(:,3);
    z(:,2) = C(:,4);
    z(:,3) = C(:,5);
    
    %% Plots over time
    %x=15:30:length(SP);
    x=1:90:length(SP);
    y=x/365;
    lstd=length(SP);
    
    %% Mean consumption level
    c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
    mclev(s,:) = nanmean(c);
    
    % Zoop overconsumption
    Zcon(s,:) = nansum(z)/lstd;
    
    %% Benthic invert biomass
    figure(11)
    subplot(4,1,1)
    plot(y,log10(C(x,1)),'b','Linewidth',1); hold on;
    plot(y,log10(C(x,2)),'k','Linewidth',1); hold on;
    plot(y,log10(C(x,5)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Benthic inverts')
    xlabel('Time (y)')
    ylabel('log10')
    legend('biomass','eaten','frac')
    legend('location','south')
    
    subplot(4,1,2)
    %plot(y,log10(C(x,1)),'b','Linewidth',1); hold on;
    plot(y,C(x,1),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Benthic inverts biomass')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    %plot(y,log10(C(x,2)),'k','Linewidth',1); hold on;
    plot(y,C(x,2),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Benthic inverts consumed')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    %plot(y,log10(C(x,5)),'r','Linewidth',1); hold on;
    plot(y,C(x,5),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0 1])
    title('Fraction of total benthic inverts consumed')
    xlabel('Time (y)')
    ylabel('Fraction')
    print('-dpng',[fpath sname lname 'oneloc_B_time.png'])
    
    %% Large Pelagic
    figure(1)
    subplot(4,1,1)
    plot(y,log10(SP(x,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MP(x,1)),'r','Linewidth',1); hold on;
    plot(y,log10(LP(x,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    legend('location','south')
    
    subplot(4,1,2)
    plot(y,log10(SP(x,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(MP(x,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(LP(x,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_P_time.png'])
    
    % Forage Fish
    figure(2)
    subplot(3,1,1)
    plot(y,log10(SF(x,1)),'b','Linewidth',1); hold on;
    plot(y,MF(x,1),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Forage Fishes ' loc])
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    legend('Immature','Adults')
    legend('location','south')
    
    subplot(3,1,2)
    plot(y,log10(SF(x,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Immature')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(3,1,3)
    plot(y,log10(MF(x,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_F_time.png'])
    
    % Demersal
    figure(3)
    subplot(4,1,1)
    plot(y,log10(SD(x,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MD(x,1)),'r','Linewidth',1); hold on;
    plot(y,log10(LD(x,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Demersal Piscivores ' loc])
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    legend('location','south')
    
    subplot(4,1,2)
    plot(y,log10(SD(x,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(MD(x,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(LD(x,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    %ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_D_time.png'])
    
    %% All size classes of all
    
    figure(5)
    plot(y,log10(SF(x,1)),'Linewidth',2); hold on;
    plot(y,log10(MF(x,1)),'Linewidth',2); hold on;
    plot(y,log10(SP(x,1)),'Linewidth',2); hold on;
    plot(y,log10(MP(x,1)),'Linewidth',2); hold on;
    plot(y,log10(LP(x,1)),'Linewidth',2); hold on;
    plot(y,log10(SD(x,1)),'Linewidth',2); hold on;
    plot(y,log10(MD(x,1)),'Linewidth',2); hold on;
    plot(y,log10(LD(x,1)),'Linewidth',2); hold on;
    legend('SF','MF','SP','MP','LP','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-10 2])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Spinup ' loc])
    print('-dpng',[fpath sname lname 'oneloc_all_sizes.png'])
    
    figure(15)
    F = SF(x,1)+MF(x,1);
    P = SP(x,1)+MP(x,1)+LP(x,1);
    D = SD(x,1)+MD(x,1)+LD(x,1);
    
    plot(y,log10(F),'r','Linewidth',2); hold on;
    plot(y,log10(P),'b','Linewidth',2); hold on;
    plot(y,log10(D),'k','Linewidth',2); hold on;
    legend('F','P','D')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    ylim([-10 2])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Spinup ' loc])
    print('-dpng',[fpath sname lname 'oneloc_all_types.png'])
    
    %% Reproduction
    rep(:,1)=MF(:,1).*MF(:,18);
    rep(:,2)=LD(:,1).*LD(:,18);
    rep(:,3)=LP(:,1).*LP(:,18);
    
    figure(6)
    subplot(4,1,1)
    plot(y,rep(x,1),'r','Linewidth',1); hold on;
    plot(y,rep(x,3),'b','Linewidth',1); hold on;
    plot(y,rep(x,2),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Reproduction ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('F','P','D')
    legend('location','north')
    
    subplot(4,1,2)
    plot(y,rep(x,1),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,rep(x,3),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Pelagic Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,rep(x,2),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_rep_time.png'])
    
    %% Maturation
    m(:,1)=log10(MF(:,19));
    m(:,2)=log10(MD(:,19));
    m(:,3)=log10(MP(:,19));
    m(:,4)=log10(LD(:,19));
    m(:,5)=log10(LP(:,19));
    
    figure(7)
    subplot(3,2,1)
    plot(y,m(x,1),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title([loc ' log10 Maturation Biomass (g m^-^2)'],'HorizontalAlignment','left')
    ylabel('Forage Fishes')
    
    subplot(3,2,3)
    plot(y,m(x,2),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('M')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,5)
    plot(y,m(x,3),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('M')
    ylabel('Pelagic Piscivores')
    xlabel('Time (y)')
    
    subplot(3,2,4)
    plot(y,m(x,4),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('L')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,6)
    plot(y,m(x,5),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('L')
    xlabel('Time (y)')
    ylabel('Pelagic Piscivores')
    print('-dpng',[fpath sname lname 'oneloc_matur_time.png'])
    
    %% Predation mortality
    %P
    figure(8)
    subplot(2,3,2)
    plot(y,log10(SP(x,17)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    
    subplot(2,3,5)
    plot(y,log10(MP(x,17)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    xlabel('Time (y)')
    
    %FF
    subplot(2,3,1)
    plot(y,log10(SF(x,17)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    ylabel('log10 Biomass eaten by predators (g m^-^2)','HorizontalAlignment','right')
    
    subplot(2,3,4)
    plot(y,log10(MF(x,17)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    %ylabel('Biomass eaten by predators (g m^-^2)')
    
    %Detritivore
    subplot(2,3,3)
    plot(y,log10(SD(x,17)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    
    subplot(2,3,6)
    plot(y,log10(MD(x,17)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    xlabel('Time (y)')
    
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_pred_sub.png'])
    
    %% All consumption in subplots
    %P
    figure(9)
    subplot(3,3,2)
    plot(y,SP(x,14),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    
    subplot(3,3,5)
    plot(y,MP(x,14),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    
    subplot(3,3,8)
    plot(y,LP(x,14),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LP')
    xlabel('Time (y)')
    
    %FF
    subplot(3,3,1)
    plot(y,SF(x,14),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    
    subplot(3,3,4)
    plot(y,MF(x,14),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass Consumed (g g^-^1 m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y,SD(x,14),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    
    subplot(3,3,6)
    plot(y,MD(x,14),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    
    subplot(3,3,9)
    plot(y,LD(x,14),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LD')
    xlabel('Time (y)')
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_consump_sub.png'])
    
    %% Recruitment
    if (y(end) > 1)
        yr=1:length(SP);
        JF=SF(yr,19);
        JD=SD(yr,19);
        JP=SP(yr,19);
        FA=MF(yr,19);
        DA=LD(yr,19);
        PA=LP(yr,19);
        
        st=1:365:length(yr);
        en=365:365:length(yr);
        SPy = NaN*ones(length(st),1);
        SFy = SPy;
        SDy = SPy;
        PAy = SPy;
        FAy = SPy;
        DAy = SPy;
        for n=1:length(st)
            SPy(n) = nansum(JP(st(n):en(n)));
            SFy(n) = nansum(JF(st(n):en(n)));
            SDy(n) = nansum(JD(st(n):en(n)));
            PAy(n) = nansum(PA(st(n):en(n)));
            FAy(n) = nansum(FA(st(n):en(n)));
            DAy(n) = nansum(DA(st(n):en(n)));
        end
        %Define recruits as new adults
        R{s,1}=FAy;
        R{s,2}=PAy;
        R{s,3}=DAy;
        
        
        figure(10)
        subplot(3,1,1)
        plot(1:length(st),log10(FAy),'r','Linewidth',2); hold on;
        xlim([1 length(st)])
        ylabel('Recruits (g m^-^2)')
        title('Forage fishes')
        
        subplot(3,1,2)
        plot(1:length(st),log10(PAy),'b','Linewidth',2); hold on;
        xlim([1 length(st)])
        ylabel({loc; 'Recruits (g m^-^2)'})
        title('Pelagic piscivores')
        
        subplot(3,1,3)
        plot(1:length(st),log10(DAy),'k','Linewidth',2); hold on;
        xlim([1 length(st)])
        ylabel('Recruits (g m^-^2)')
        title('Demersal piscivores')
        print('-dpng',[fpath sname lname 'oneloc_recruitment.png'])
    end
    
end

save([dpath sname 'consump.mat'],'mclev','Zcon');
%     csvwrite([dpath sname 'clevel.csv'],mclev);
%     csvwrite([dpath sname 'Zconsump.csv'],Zcon);

if (y(end) > 1)
%     csvwrite([dpath sname 'Recruits.csv'],R);
    save([dpath sname 'Recruits.mat'],'R');
end


