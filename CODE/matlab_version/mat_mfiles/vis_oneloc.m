% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefOrig/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefOrig/';

spots = {'GB','EBS','OSP','HOT','BATS','NS'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Oneloc_hist_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),2);
R=NaN*ones(30,3,length(spots));
for s=1:length(spots)
    %%
    close all
    loc = spots{s};
    %lname = [loc '_'];
    lname = ['phen_' loc '_'];
    
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    C = csvread([dpath sname lname 'Cobalt.csv']);
    z(:,1) = C(:,2);
    z(:,2) = C(:,3);
    
    %% Plots over time
    x=1:length(SP);
    y=x/365;
    lstd=length(SP);
    id1 = 0:365:(lstd-1);
    id2 = 365:365:(lstd);
    ID  = [id1 id2];
    
    %% Mean consumption level
    c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
    mclev(s,:) = nanmean(c);
    
    % Zoop overconsumption
    Zcon(s,:) = nansum(z)/lstd;
    
    
    %% Piscivore
    figure(1)
    subplot(4,1,1)
    plot(y((end-(30*365)):end),SP((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    plot(y((end-(30*365)):end),MP((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    plot(y((end-(30*365)):end),LP((end-(30*365)):end,1),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title(['Historical Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y((end-(30*365)):end),SP((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y((end-(30*365)):end),MP((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y((end-(30*365)):end),LP((end-(30*365)):end,1),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_pisc_time.png'])
    
    %% Planktivore
    figure(2)
    subplot(3,1,1)
    plot(y((end-(30*365)):end),SF((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    plot(y((end-(30*365)):end),MF((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title(['Historical Forage Fishes ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('Immature','Adults')
    
    subplot(3,1,2)
    plot(y((end-(30*365)):end),SF((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Immature')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,1,3)
    plot(y((end-(30*365)):end),MF((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    print('-dpng',[fpath sname lname 'oneloc_plan_time.png'])
    
    %% Detritivore
    figure(3)
    subplot(4,1,1)
    plot(y((end-(30*365)):end),SD((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    plot(y((end-(30*365)):end),MD((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    plot(y((end-(30*365)):end),LD((end-(30*365)):end,1),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title(['Historical Demersal Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y((end-(30*365)):end),SD((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y((end-(30*365)):end),MD((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y((end-(30*365)):end),LD((end-(30*365)):end,1),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_detr_time.png'])
    
    %% All biomass in subplots
    %SP
    figure(4)
    subplot(3,3,2)
    plot(y((end-(30*365)):end),SP((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title({loc; 'SP'})
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,5)
    plot(y((end-(30*365)):end),MP((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MP')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,8)
    plot(y((end-(30*365)):end),LP((end-(30*365)):end,1),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('LP')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    %FF
    subplot(3,3,1)
    plot(y((end-(30*365)):end),SF((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('SF')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,4)
    plot(y((end-(30*365)):end),MF((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y((end-(30*365)):end),SD((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('SD')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,6)
    plot(y((end-(30*365)):end),MD((end-(30*365)):end,1),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MD')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,9)
    plot(y((end-(30*365)):end),LD((end-(30*365)):end,1),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('LD')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_sub.png'])
    
    %% All size classes of all
    
    figure(5)
    plot(y((end-(30*365)):end),SP((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),MP((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),LP((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),SF((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),MF((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),SD((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),MD((end-(30*365)):end,1),'Linewidth',1); hold on;
    plot(y((end-(30*365)):end),LD((end-(30*365)):end,1),'Linewidth',1); hold on;
    legend('SP','MP','LP','SF','MF','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(end-(30*365)) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    title(['Historical ' loc])
    print('-dpng',[fpath sname lname 'oneloc_all_sizes.png'])
    
    %% Reproduction
    rep(:,1)=MF(:,1).*MF(:,18);
    rep(:,2)=LD(:,1).*LD(:,18);
    rep(:,3)=LP(:,1).*LP(:,18);
    
    figure(6)
    subplot(4,1,1)
    plot(y((end-(30*365)):end),rep((end-(30*365)):end,:),'Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title(['Historical Reproduction ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('F','D','SP')
    
    subplot(4,1,2)
    plot(y((end-(30*365)):end),rep((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y((end-(30*365)):end),rep((end-(30*365)):end,2),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y((end-(30*365)):end),rep((end-(30*365)):end,3),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('Pelagic Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_rep_time.png'])
    
    %% Maturation
    m(:,1)=MF(:,19);
    m(:,2)=MD(:,19);
    m(:,3)=MP(:,19);
    m(:,4)=LD(:,19);
    m(:,5)=LP(:,19);
    
    figure(7)
    subplot(3,2,1)
    plot(y((end-(30*365)):end),m((end-(30*365)):end,1),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title([loc ' Maturation Biomass (g m^-^2)'],'HorizontalAlignment','left')
    ylabel('Forage Fishes')
    
    subplot(3,2,3)
    plot(y((end-(30*365)):end),m((end-(30*365)):end,2),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('M')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,5)
    plot(y((end-(30*365)):end),m((end-(30*365)):end,3),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('M')
    ylabel('Pelagic Piscivores')
    xlabel('Time (y)')
    
    subplot(3,2,4)
    plot(y((end-(30*365)):end),m((end-(30*365)):end,4),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('L')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,6)
    plot(y((end-(30*365)):end),m((end-(30*365)):end,5),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('L')
    xlabel('Time (y)')
    ylabel('Pelagic Piscivores')
    print('-dpng',[fpath sname lname 'oneloc_matur_time.png'])
    
    %% Predation mortality
    %SP
    figure(8)
    subplot(2,3,2)
    plot(y((end-(30*365)):end),SP((end-(30*365)):end,17),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title({loc; 'SP'})
    
    subplot(2,3,5)
    plot(y((end-(30*365)):end),MP((end-(30*365)):end,17),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MP')
    xlabel('Time (y)')
    
    %FF
    subplot(2,3,1)
    plot(y((end-(30*365)):end),SF((end-(30*365)):end,17),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('SF')
    ylabel('Biomass eaten by predators (g m^-^2)','HorizontalAlignment','right')
    
    subplot(2,3,4)
    plot(y((end-(30*365)):end),MF((end-(30*365)):end,17),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MF')
    xlabel('Time (y)')
    %ylabel('Biomass eaten by predators (g m^-^2)')
    
    %Detritivore
    subplot(2,3,3)
    plot(y((end-(30*365)):end),SD((end-(30*365)):end,17),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('SD')
    
    subplot(2,3,6)
    plot(y((end-(30*365)):end),MD((end-(30*365)):end,17),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MD')
    xlabel('Time (y)')
    
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_pred_sub.png'])
    
    %% All consumption in subplots
    %SP
    figure(9)
    subplot(3,3,2)
    plot(y((end-(30*365)):end),SP((end-(30*365)):end,14),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title({loc; 'SP'})
    
    subplot(3,3,5)
    plot(y((end-(30*365)):end),MP((end-(30*365)):end,14),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MP')
    
    subplot(3,3,8)
    plot(y((end-(30*365)):end),LP((end-(30*365)):end,14),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('LP')
    xlabel('Time (y)')
    
    %FF
    subplot(3,3,1)
    plot(y((end-(30*365)):end),SF((end-(30*365)):end,14),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('SF')
    
    subplot(3,3,4)
    plot(y((end-(30*365)):end),MF((end-(30*365)):end,14),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass Consumed (g m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y((end-(30*365)):end),SD((end-(30*365)):end,14),'b','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('SD')
    
    subplot(3,3,6)
    plot(y((end-(30*365)):end),MD((end-(30*365)):end,14),'r','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('MD')
    
    subplot(3,3,9)
    plot(y((end-(30*365)):end),LD((end-(30*365)):end,14),'k','Linewidth',1); hold on;
    xlim([y(end-(30*365)) y(end)])
    title('LD')
    xlabel('Time (y)')
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_consump_sub.png'])
    
    %% Recruitment
    t=1:length(y);
    yr=t((end-(30*365)+1):end);
    SFL=m(yr,3);
    SDL=m(yr,4);
    SPL=m(yr,5);
    FA=MF(yr,1);
    DA=LD(yr,1);
    PA=LP(yr,1);
    
    st=1:365:length(yr);
    en=365:365:length(yr);
    SPy = NaN*ones(30,1);
    SFy = SPy;
    SDy = SPy;
    PAy = SPy;
    FAy = SPy;
    DAy = SPy;
    for n=1:30
        SPy(n) = nansum(SPL(st(n):en(n)));
        SFy(n) = nansum(SFL(st(n):en(n)));
        SDy(n) = nansum(SDL(st(n):en(n)));
        PAy(n) = nansum(PA(st(n):en(n)));
        FAy(n) = nansum(FA(st(n):en(n)));
        DAy(n) = nansum(DA(st(n):en(n)));
    end
    
    R(:,1,s)=SFy;
    R(:,2,s)=SPy;
    R(:,3,s)=SDy;
    
    %
    figure(10)
    subplot(3,1,1)
    plot(1976:2005,SPy,'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel({loc; 'Recruits (g km^-^2)'})
    title('Pelagic piscivores')
    
    subplot(3,1,2)
    plot(1976:2005,SFy,'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Recruits (g km^-^2)')
    title('Forage fishes')
    
    subplot(3,1,3)
    plot(1976:2005,SDy,'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Recruits (g km^-^2)')
    title('Demersal piscivores')
    print('-dpng',[fpath sname lname 'oneloc_recruitment.png'])
end

if (length(lname) > 5)
    save([dpath sname lname(1:5) 'consump.mat'],'mclev','Zcon');
    csvwrite([dpath sname lname(1:5) 'clevel.csv'],mclev);
    csvwrite([dpath sname lname(1:5) 'Zconsump.csv'],Zcon);
    csvwrite([dpath sname lname(1:5) 'Recruits.csv'],R);
    save([dpath sname lname(1:5) 'Recruits.mat'],'R');
else
    save([dpath sname 'consump.mat'],'mclev','Zcon');
    csvwrite([dpath sname 'clevel.csv'],mclev);
    csvwrite([dpath sname 'Zconsump.csv'],Zcon);
    csvwrite([dpath sname 'Recruits.csv'],R);
    save([dpath sname 'Recruits.mat'],'R');
end


