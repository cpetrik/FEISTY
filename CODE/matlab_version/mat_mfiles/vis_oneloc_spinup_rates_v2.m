% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 50 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
sreff = {'1.0','0.5','0.1','0.05','0.01'};
encs = logspace(2,4,10); %linspace(10,90,9); %
cmaxs = linspace(10,100,10);
efn = 70;
tefn = num2str(efn);
% cfn = 20;
% tcfn = num2str(cfn);
fcrit = 20;
nmort = '1';
kad = 50;
D = 'D050';
J = 'J050';
Sm = 'Sm025';
BE = '05';
CC = '050';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
sname2 = '';

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

mass = [M_s;M_m;M_l];
mass = repmat(mass,1,length(spots));
L = [L_s;L_m;L_l];

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%%
% for E = 1:length(encs)
%     efn = encs(E);
%     tefn = num2str(round(efn));
    for C = 1:length(cmaxs)
        cfn = cmaxs(C);
        tcfn = num2str(cfn);
    for R = 3;%1:length(RE)
        rfrac = RE{R};
        close all
        dp = ['Dc_enc',tefn,'_cmax-metab',tcfn,'_fcrit',num2str(fcrit),'_',D,'_',J,'_',Sm,...
            '_nmort',nmort,'_BE',BE,'_CC',CC,'_RE',rfrac];
        dpath = [datap char(dp) '/'];
        fpath = [figp char(dp) '/'];
        if (~isdir([figp char(dp)]))
            mkdir([figp char(dp)])
        end
        cfile = char(dp);
        
        Psum = NaN*ones(3,length(spots));
        Fsum = NaN*ones(2,length(spots));
        Dsum = Psum;
        Bmean = NaN*ones(1,length(spots));
        Pmean = Psum;
        Fmean = Fsum;
        Dmean = Psum;
        Pmgr = Psum;
        Fmgr = Fsum;
        Dmgr = Psum;
        Pgge = Psum;
        Fgge = Fsum;
        Dgge = Psum;
        Pprod = Psum;
        Fprod = Fsum;
        Dprod = Psum;
        Pcon = Psum;
        Fcon = Fsum;
        Dcon = Psum;
        Prep = Fsum;
        Frep = Fsum;
        Drep = Fsum;
        Pmet = Psum;
        Fmet = Fsum;
        Dmet = Psum;
        Ppred = Psum;
        Fpred = Fsum;
        Dpred = Psum;
        Pnat = Psum;
        Fnat = Fsum;
        Dnat = Psum;
        Pfish = Psum;
        Ffish = Fsum;
        Dfish = Psum;
        Ptotcatch = Psum;
        Ftotcatch = Fsum;
        Dtotcatch = Psum;
        Plev = Psum;
        Flev = Fsum;
        Dlev = Psum;
        all_mean=NaN*ones(3,4,length(spots));
        z = NaN*ones(length(spots),3);
        
        %%
        for s=1:length(spots)
            %%
            loc = spots{s};
            lname = [loc '_'];
            SP = csvread([dpath sname lname 'Sml_p.csv']);
            SF = csvread([dpath sname lname 'Sml_f.csv']);
            SD = csvread([dpath sname lname 'Sml_d.csv']);
            MP = csvread([dpath sname lname 'Med_p.csv']);
            MF = csvread([dpath sname lname 'Med_f.csv']);
            MD = csvread([dpath sname lname 'Med_d.csv']);
            LP = csvread([dpath sname lname 'Lrg_p.csv']);
            LD = csvread([dpath sname lname 'Lrg_d.csv']);
            CO = csvread([dpath sname lname 'Cobalt.csv']);
            
            t=1:length(SP);
            lyr=t((end-12+1):end);
            
            %% Final mean biomass in each size
            
            SP_mean=mean(SP(lyr,1));
            SF_mean=mean(SF(lyr,1));
            SD_mean=mean(SD(lyr,1));
            MP_mean=mean(MP(lyr,1));
            MF_mean=mean(MF(lyr,1));
            MD_mean=mean(MD(lyr,1));
            LP_mean=mean(LP(lyr,1));
            LD_mean=mean(LD(lyr,1));
            B_mean =mean(CO(lyr,1));
            
            P_mean=[SP_mean;MP_mean;LP_mean];
            F_mean=[SF_mean;MF_mean];
            D_mean=[SD_mean;MD_mean;LD_mean];
            
            Pmean(:,s) = P_mean;
            Fmean(:,s) = F_mean;
            Dmean(:,s) = D_mean;
            Bmean(:,s) = B_mean;
            
            all_mean(1:2,1,s) = F_mean;
            all_mean(:,2,s) = P_mean;
            all_mean(:,3,s) = D_mean;
            all_mean(1,4,s) = B_mean;
            
            %% Growth rate (nu - energy for biomass production)
            SP_mgr=nanmean(SP(lyr,15));
            SF_mgr=nanmean(SF(lyr,15));
            SD_mgr=nanmean(SD(lyr,15));
            MP_mgr=nanmean(MP(lyr,15));
            MF_mgr=nanmean(MF(lyr,15));
            MD_mgr=nanmean(MD(lyr,15));
            LP_mgr=nanmean(LP(lyr,15));
            LD_mgr=nanmean(LD(lyr,15));
            
            P_mgr=[SP_mgr;MP_mgr;LP_mgr];
            F_mgr=[SF_mgr;MF_mgr];
            D_mgr=[SD_mgr;MD_mgr;LD_mgr];
            
            Pmgr(:,s) = P_mgr;
            Fmgr(:,s) = F_mgr;
            Dmgr(:,s) = D_mgr;
            
            %% Consump per biomass (I)
            SP_con=nanmean(SP(lyr,14));
            SF_con=nanmean(SF(lyr,14));
            SD_con=nanmean(SD(lyr,14));
            MP_con=nanmean(MP(lyr,14));
            MF_con=nanmean(MF(lyr,14));
            MD_con=nanmean(MD(lyr,14));
            LP_con=nanmean(LP(lyr,14));
            LD_con=nanmean(LD(lyr,14));
            
            P_con=[SP_con;MP_con;LP_con];
            F_con=[SF_con;MF_con];
            D_con=[SD_con;MD_con;LD_con];
            
            Pcon(:,s) = P_con;
            Fcon(:,s) = F_con;
            Dcon(:,s) = D_con;
            
            %% Feeding level = con/cmax
            SP_lev=nanmean(SP(lyr,20));
            SF_lev=nanmean(SF(lyr,20));
            SD_lev=nanmean(SD(lyr,20));
            MP_lev=nanmean(MP(lyr,20));
            MF_lev=nanmean(MF(lyr,20));
            MD_lev=nanmean(MD(lyr,20));
            LP_lev=nanmean(LP(lyr,20));
            LD_lev=nanmean(LD(lyr,20));
            
            P_lev=[SP_lev;MP_lev;LP_lev];
            F_lev=[SF_lev;MF_lev];
            D_lev=[SD_lev;MD_lev;LD_lev];
            
            Plev(:,s) = P_lev;
            Flev(:,s) = F_lev;
            Dlev(:,s) = D_lev;
            
            %% Fraction zoop losses consumed
            z(s,1) = nanmean(CO(lyr,3));
            z(s,2) = nanmean(CO(lyr,4));
            z(s,3) = nanmean(CO(lyr,5));
            
            %% Size spectrum (sum stages)
            spec = nansum(all_mean(:,:,s),2);
            
            %% Production (= nu * biom)
            SP_prod=nanmean(SP(lyr,21));
            SF_prod=nanmean(SF(lyr,21));
            SD_prod=nanmean(SD(lyr,21));
            MP_prod=nanmean(MP(lyr,21));
            MF_prod=nanmean(MF(lyr,21));
            MD_prod=nanmean(MD(lyr,21));
            LP_prod=nanmean(LP(lyr,21));
            LD_prod=nanmean(LD(lyr,21));
            
            P_prod=[SP_prod;MP_prod;LP_prod];
            F_prod=[SF_prod;MF_prod];
            D_prod=[SD_prod;MD_prod;LD_prod];
            
            Pprod(:,s) = P_prod;
            Fprod(:,s) = F_prod;
            Dprod(:,s) = D_prod;
            
            %% Reproduction
            F_rep(1,1)=nanmean(MF(lyr,18));
            D_rep(1,1)=nanmean(LD(lyr,18));
            P_rep(1,1)=nanmean(LP(lyr,18));
            F_rep(2,1)=nanmean(MF(lyr,1).*MF(lyr,18));
            D_rep(2,1)=nanmean(LD(lyr,1).*LD(lyr,18));
            P_rep(2,1)=nanmean(LP(lyr,1).*LP(lyr,18));
            
            Prep(:,s) = P_rep;
            Frep(:,s) = F_rep;
            Drep(:,s) = D_rep;
            
            %% Metabolism
            SP_met=nanmean(SP(lyr,24));
            SF_met=nanmean(SF(lyr,24));
            SD_met=nanmean(SD(lyr,24));
            MP_met=nanmean(MP(lyr,24));
            MF_met=nanmean(MF(lyr,24));
            MD_met=nanmean(MD(lyr,24));
            LP_met=nanmean(LP(lyr,24));
            LD_met=nanmean(LD(lyr,24));
            
            P_met=[SP_met;MP_met;LP_met];
            F_met=[SF_met;MF_met];
            D_met=[SD_met;MD_met;LD_met];
            
            Pmet(:,s) = P_met;
            Fmet(:,s) = F_met;
            Dmet(:,s) = D_met;
            
            %% Predation
            SP_pred=nanmean(SP(lyr,22));
            SF_pred=nanmean(SF(lyr,22));
            SD_pred=nanmean(SD(lyr,22));
            MP_pred=nanmean(MP(lyr,22));
            MF_pred=nanmean(MF(lyr,22));
            MD_pred=nanmean(MD(lyr,22));
            LP_pred=nanmean(LP(lyr,22));
            LD_pred=nanmean(LD(lyr,22));
            
            P_pred=[SP_pred;MP_pred;LP_pred];
            F_pred=[SF_pred;MF_pred];
            D_pred=[SD_pred;MD_pred;LD_pred];
            
            Ppred(:,s) = P_pred;
            Fpred(:,s) = F_pred;
            Dpred(:,s) = D_pred;
            
            %% Natural mortality
            Pnat(1,s)=nanmean(SP(lyr,23));
            Fnat(1,s)=nanmean(SF(lyr,23));
            Dnat(1,s)=nanmean(SD(lyr,23));
            Pnat(2,s)=nanmean(MP(lyr,23));
            Fnat(2,s)=nanmean(MF(lyr,23));
            Dnat(2,s)=nanmean(MD(lyr,23));
            Pnat(3,s)=nanmean(LP(lyr,23));
            Dnat(3,s)=nanmean(LD(lyr,23));
            
            %% Fishing
            MP_fish=nanmean(MP(lyr,25));
            MF_fish=nanmean(MF(lyr,25));
            MD_fish=nanmean(MD(lyr,25));
            LP_fish=nanmean(LP(lyr,25));
            LD_fish=nanmean(LD(lyr,25));
            
            P_fish=[0;MP_fish;LP_fish];
            F_fish=[0;MF_fish];
            D_fish=[0;MD_fish;LD_fish];
            
            Pfish(:,s) = P_fish;
            Ffish(:,s) = F_fish;
            Dfish(:,s) = D_fish;
            
            MP_totcatch=nansum(MP(lyr,25));
            MF_totcatch=nansum(MF(lyr,25));
            MD_totcatch=nansum(MD(lyr,25));
            LP_totcatch=nansum(LP(lyr,25));
            LD_totcatch=nansum(LD(lyr,25));
            
            P_totcatch=[0;MP_totcatch;LP_totcatch];
            F_totcatch=[0;MF_totcatch];
            D_totcatch=[0;MD_totcatch;LD_totcatch];
            
            Ptotcatch(:,s) = P_totcatch;
            Ftotcatch(:,s) = F_totcatch;
            Dtotcatch(:,s) = D_totcatch;
            
            %% Total mortality w/o fishing
            %model lengths
            L(1) = 10;
            L(2) = 200;
            L(3) = 1e3;
            %model mass in grams
            M = 0.01 .* (0.1.*L).^3;
            %Andersen & Beyer mortality rate per year (natural + predation)
            %physiol mort * growth constant * M^-0.25
            AB = (0.35 .* 4.5 .* M.^(-0.25)) ./365;
            
            Fmort = Fpred + Fnat;
            Pmort = Ppred + Pnat;
            Dmort = Dpred + Dnat;
            
            %% Total mortality w/ fishing
            Fmortf = Fpred + Fnat + Ffish;
            Pmortf = Ppred + Pnat + Pfish;
            Dmortf = Dpred + Dnat + Dfish;
            
            %% Gross growth efficiency (= nu/consump)
            SP_gge=nanmean(SP(lyr,15)./SP(lyr,14));
            SF_gge=nanmean(SF(lyr,15)./SF(lyr,14));
            SD_gge=nanmean(SD(lyr,15)./SD(lyr,14));
            MP_gge=nanmean(MP(lyr,15)./MP(lyr,14));
            MF_gge=nanmean(MF(lyr,15)./MF(lyr,14));
            MD_gge=nanmean(MD(lyr,15)./MD(lyr,14));
            LP_gge=nanmean(LP(lyr,15)./LP(lyr,14));
            LD_gge=nanmean(LD(lyr,15)./LD(lyr,14));
            
            P_gge=[SP_gge;MP_gge;LP_gge];
            F_gge=[SF_gge;MF_gge];
            D_gge=[SD_gge;MD_gge;LD_gge];
            
            Pgge(:,s) = P_gge;
            Fgge(:,s) = F_gge;
            Dgge(:,s) = D_gge;
            
        end %spots
        save([dpath sname 'lastyr_sum_mean_biom.mat'],'Psum','Fsum',...
            'Dsum','Pmean','Fmean','Dmean','all_mean',...
            'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
            'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
            'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
            'Dtotcatch','Pgge','Fgge','Dgge','Plev','Flev','Dlev','Bmean');
        
        %% PLOTS
        mlev = [Flev;Plev;Dlev];
        for s=1:length(spots)
            loc = spots{s};
            lname = [loc '_'];
            
            %% Biomass
            f21 = figure(21);
            subplot(4,3,s)
            plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 6])
            ylim([-3 2])
            set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
            if (s==4)
                ylabel('log10 Mean Biom (g m^-^2) in final year')
            end
            title(loc)
            if (s==3)
                stamp(cfile)
            end
            
            %% Feeding level
            f2=figure(2);
            subplot(4,3,s)
            bar(mlev(:,s),'k')
            ylim([0 1])
            xlim([0 9])
            set(gca,'XTickLabel',[]);
            for n=1:8
                text(n-0.5,-0.2,stages{n},'Rotation',45)
            end
            title(spots{s})
            if (s==4)
                ylabel('Feeding level')
            end
            if (s==11)
                subplot(4,3,12)
                bar(mean(mlev,2),'k')
                ylim([0 1])
                xlim([0 9])
                set(gca,'XTickLabel',[]);
                for n=1:8
                    text(n-0.5,-0.2,stages{n},'Rotation',45)
                end
                title('Mean all')
                stamp(cfile)
            end
            
            %% Growth rate (nu - energy for biomass production)
            f3 = figure(3);
            subplot(2,2,1)
            plot(s-0.25,Fmgr(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Pmgr(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dmgr(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean growth rate (g g^-^1 d^-^1) in final year')
            title('S')
            
            subplot(2,2,2)
            plot(s-0.25,(Fmgr(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Pmgr(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmgr(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean growth/repro rate (g g^-^1 d^-^1) in final year')
            title('M')
            
            subplot(2,2,3)
            plot(s,(Pmgr(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmgr(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
            title('L')
            if (s==3)
                stamp(cfile)
            end
            
            %% Fraction zoop losses consumed
            f5 = figure(5);
            subplot(4,3,s)
            bar(z(s,:),'k'); hold on;
            xlim([0 4])
            ylim([0 1])
            set(gca,'XTick',1:3,'XTickLabel',{'MZ','LZ','Bent'})
            if (s==4)
                ylabel('Fraction flux consumed')
            end
            title(loc)
            if (s==3)
                stamp(cfile)
            end
            
            %% Production (= nu * biom)
            f8 = figure(8);
            subplot(2,2,1)
            plot(s-0.25,Fprod(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Pprod(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dprod(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
            title('S')
            
            subplot(2,2,2)
            plot(s-0.25,(Fprod(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Pprod(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dprod(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
            title('M')
            
            subplot(2,2,3)
            plot(s,(Pprod(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dprod(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
            title('L')
            if (s==1)
                stamp(cfile)
            end
            
            %% Reproduction
            f9 = figure(9);
            subplot(1,2,1)
            plot(s-0.25,Frep(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Prep(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Drep(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.01,spots{n},'Rotation',45)
            end
            ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
            
            subplot(1,2,2)
            plot(s-0.25,(Frep(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Prep(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Drep(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.01,spots{n},'Rotation',45)
            end
            ylabel('Mean biom reproduced (g d^-^1) in final year')
            if (s==1)
                stamp(cfile)
            end
            
            %% Metabolism
            f10 = figure(10);
            subplot(2,2,1)
            plot(s-0.25,Fmet(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Pmet(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dmet(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
            title('S')
            
            subplot(2,2,2)
            plot(s-0.25,(Fmet(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Pmet(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmet(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
            title('M')
            
            subplot(2,2,3)
            plot(s,(Pmet(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmet(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
            title('L')
            if (s==3)
                stamp(cfile)
            end
            
            %% Predation
            f11 = figure(11);
            subplot(1,2,1)
            plot(s-0.25,Fpred(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Ppred(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dpred(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.001,spots{n},'Rotation',45)
            end
            ylabel('Mean predation rate (d^-^1) in final year')
            title('S')
            
            subplot(1,2,2)
            plot(s-0.25,(Fpred(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Ppred(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dpred(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.001,spots{n},'Rotation',45)
            end
            ylabel('Mean predation rate (d^-^1) in final year')
            title('M')
            if (s==3)
                stamp(cfile)
            end
            
            %% Fishing
            f12 = figure(12);
            subplot(1,2,1)
            plot(s-0.25,Ftotcatch(2,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Ptotcatch(2,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dtotcatch(2,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.1,spots{n},'Rotation',45)
            end
            ylabel('Total catch (g) in final year')
            
            subplot(1,2,2)
            plot(s,Ptotcatch(3,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dtotcatch(3,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.1,spots{n},'Rotation',45)
            end
            ylabel('Total catch (g) in final year')
            if (s==1)
                stamp(cfile)
            end
            
            %% Total mortality w/o fishing
            %model mass in grams
            M = 0.01 .* (0.1.*L).^3;
            %Andersen & Beyer mortality rate per year (natural + predation)
            %physiol mort * growth constant * M^-0.25
            AB = (0.35 .* 4.5 .* M.^(-0.25)) ./365;
            
            Fmort = Fpred + Fnat;
            Pmort = Ppred + Pnat;
            Dmort = Dpred + Dnat;
            
            f13=figure(13);
            subplot(2,2,1)
            plot(s-0.25,Fmort(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Pmort(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dmort(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            if(s==11)
                plot(0:10,AB(1)*ones(11,1),'--k'); hold on;
            end
            ylabel('Mean mortality rate w/o fishing (d^-^1) in final year')
            title('S')
            
            subplot(2,2,2)
            plot(s-0.25,(Fmort(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Pmort(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmort(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            if(s==11)
                plot(0:10,AB(2)*ones(11,1),'--k'); hold on;
            end
            title('M')
            
            subplot(2,2,3)
            plot(s,(Pmort(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmort(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            if(s==11)
                plot(0:10,AB(3)*ones(11,1),'--k'); hold on;
            end
            title('L')
            
            %% Total mortality w/ fishing
            Fmortf = Fpred + Fnat + Ffish;
            Pmortf = Ppred + Pnat + Pfish;
            Dmortf = Dpred + Dnat + Dfish;
            
            f14=figure(14);
            subplot(2,2,1)
            plot(s-0.25,Fmortf(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Pmortf(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dmortf(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            if(s==11)
                plot(0:10,AB(1)*ones(11,1),'--k'); hold on;
            end
            ylabel('Mean mortality rate w/fishing (d^-^1) in final year')
            title('S')
            
            subplot(2,2,2)
            plot(s-0.25,(Fmortf(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Pmortf(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmortf(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            if(s==11)
                plot(0:10,AB(2)*ones(11,1),'--k'); hold on;
            end
            title('M')
            
            subplot(2,2,3)
            plot(s,(Pmortf(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dmortf(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,0,spots{n},'Rotation',45)
            end
            if(s==11)
                plot(0:10,AB(3)*ones(11,1),'--k'); hold on;
            end
            title('L')
            
            %% Gross growth efficiency (= nu/consump)
            f15 = figure(15);
            subplot(2,2,1)
            plot(s-0.25,Fgge(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,Pgge(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,Dgge(1,s),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            ylim([-0.5 0.5])
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.6,spots{n},'Rotation',45)
            end
            ylabel('Mean gross growth efficiency in final year')
            title('S')
            
            subplot(2,2,2)
            plot(s-0.25,(Fgge(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(s,(Pgge(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dgge(2,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            ylim([-0.5 0.5])
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.6,spots{n},'Rotation',45)
            end
            ylabel('Mean gross growth efficiency in final year')
            title('M')
            
            subplot(2,2,3)
            plot(s,(Pgge(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(s+0.25,(Dgge(3,s)),'sk',...
                'MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            ylim([-0.5 0.5])
            xlim([0 12])
            set(gca,'XTickLabel',[]);
            for n=1:11
                text(n-0.5,-0.6,spots{n},'Rotation',45)
            end
            ylabel('Mean gross growth efficiency in final year')
            title('L')
            if (s==3)
                stamp(cfile)
            end
            
        end
        print(f21,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass_axes.png'])
        print(f2,'-dpng',[fpath sname sname2 'All_oneloc_con_level.png'])
        print(f3,'-dpng',[fpath sname sname2 'All_oneloc_nu.png'])
        print(f5,'-dpng',[fpath sname sname2 'All_oneloc_frac_zoop_loss.png'])
        print(f8,'-dpng',[fpath sname sname2 'All_oneloc_prod.png'])
        print(f9,'-dpng',[fpath sname sname2 'All_oneloc_rep.png'])
        print(f10,'-dpng',[fpath sname sname2 'All_oneloc_met.png'])
        print(f11,'-dpng',[fpath sname sname2 'All_oneloc_pred.png'])
        print(f12,'-dpng',[fpath sname sname2 'All_oneloc_catch.png'])
        print(f13,'-dpng',[fpath sname sname2 'All_oneloc_mort_nof.png'])
        print(f14,'-dpng',[fpath sname sname2 'All_oneloc_mort_f.png'])
        print(f15,'-dpng',[fpath sname sname2 'All_oneloc_gge.png'])
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        figure(16);
        plot((1-0.2):11,log10(fishsp(1,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:11,log10(fishsp(2,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot((1+0.2):12,log10(fishsp(3,:)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 12])
        ylim([-2 2])
        set(gca,'XTick',1:11,'XTickLabel',[])
        for n=1:11
            text(n,-2.2,spots{n},'HorizontalAlignment','center')
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title('All stages')
        stamp(cfile)
        print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_type.png'])
        
        sumspec = squeeze(nansum(nansum(all_mean)));
        
        figure(17);
        plot(1:11,log10(sumspec),'k.','MarkerSize',25); hold on;
        xlim([0 12])
        ylim([-2 2])
        set(gca,'XTick',1:11,'XTickLabel',[])
        for n=1:11
            text(n,-2.1,spots{n},'HorizontalAlignment','center')
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title('All fishes and stages')
        stamp(cfile)
        print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec.png'])
        
        %% Consump g/g/d --> g/d --> g/y
        Pcon = Pcon .* mass .* 365;
        Fcon = Fcon .* mass(1:2,:) .* 365;
        Dcon = Dcon .* mass .* 365;
        
        f18 = figure(18);
        subplot(3,1,1)
        plot(0.75:1:10.75,Fcon(1,:),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:1:11,Pcon(1,:),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1.25:1:11.25,Dcon(1,:),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 12])
        ylabel('S')
        set(gca,'XTick',1:11,'XTickLabel',spots)
        title('Mean consumption (g y^-^1) in final year')
        
        subplot(3,1,2)
        plot(0.75:1:10.75,Fcon(2,:),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:1:11,Pcon(2,:),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1.25:1:11.25,Dcon(2,:),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 12])
        ylabel('M')
        set(gca,'XTick',1:11,'XTickLabel',spots)
        
        subplot(3,1,3)
        plot(1:1:11,Pcon(3,:),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1.25:1:11.25,Dcon(3,:),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 12])
        ylabel('L')
        set(gca,'XTick',1:11,'XTickLabel',spots)
        xlabel('Location')
        stamp(cfile)
        
        print(f18,'-dpng',[fpath sname sname2 'All_oneloc_consump_gyr.png'])
        
        %% Consump vs. weight
        A = 4.39;
        fc = 0.2;
        f0 = 0.6;
        epsassim = 0.6;
        n = 3/4;
        
        w = logspace(-3, 5);
        
        AvailEnergy = A*w.^n;
        Consumption = A / (epsassim*(f0-fc)) * w.^n;
        
        f19=figure(19);
        for s=1:length(spots)
            subplot(1,3,1)
            loglog((mass(1:2,1)),(Fcon(:,s)),'.',...
                'Color',cm{s},'MarkerSize',25); hold on;
            title('F')
            xlabel('Mass (g)')
            ylabel('Mean consumption (g y^-^1) in final year')
            legend(spots)
            legend('location','northwest')
            %axis([-5 5 -1 5])
            
            subplot(1,3,2)
            loglog((mass(:,1)),(Pcon(:,s)),'.',...
                'Color',cm{s},'MarkerSize',25); hold on;
            title('P')
            xlabel('Mass (g)')
            %axis([-5 5 -1 5])
            
            subplot(1,3,3)
            loglog((mass(:,1)),(Dcon(:,s)),'.',...
                'Color',cm{s},'MarkerSize',25); hold on;
            title('D')
            xlabel('Mass (g)')
            %axis([-5 5 -1 5])
            stamp(cfile)
        end
        subplot(1,3,1)
        loglog(w, Consumption,'k')
        
        subplot(1,3,2)
        loglog(w, Consumption,'k')
        
        subplot(1,3,3)
        loglog(w, Consumption,'k')
        print(f19,'-dpng',[fpath sname 'All_oneloc_consump_gyr_vs_weight_compare.png'])
        
    end %RE
end %cmax
% end %enc
