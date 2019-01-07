% Visualize output of POEM mortality tests
% Spinup at one location
% 50 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_big_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Fishing/';

fnum = [0.1:0.1:1];%,1.2:0.2:2];
fsim = {'.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};%,'1.2','1.4','1.6','1.8','2'};
frate = {'01','02','03','04','05','06','07','08','09','10'};%,'12','14','16','18','20'};
RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
% nmrt = [0,2:5];
% Mort = {'None','Hartvig','mizer','J&C','P&W'};
nmort = '3';
fcrit = 40;
kad = 100;
pref = 'D100';
BE = '05';
CC = '050';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

%%
for r = 1:length(RE)
    close all
    rfrac = RE{r};
    cfile2 = ['Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit' num2str(fcrit) ...
        '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac ...
        '_LD_fishing_catch'];
    for i=1:length(frate)
        F = frate{i};
        dp = ['Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LD_fish' F];
        dpath = [datap char(dp) '/'];
        %fpath = [figp char(dp) '/'];
        fpath = figp;
        cfile = char(dp);

        
        %%
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
            
        end
        
        save([dpath sname 'lastyr_sum_mean_biom'],'Pmean','Fmean','Dmean','all_mean',...
            'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
            'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
            'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
            'Dtotcatch','Pgge','Fgge','Dgge');
        
    end
    
    %%
%     for s=1:length(spots)
%         
%         loc = spots{s};
%         lname = [loc '_'];
%         close all
%         %%
%         for i=1:length(frate)
%             F = frate{i};
%             dp = ['Dc_TrefO_cmax-metab2_enc1_MFeqMP_fcrit' num2str(fcrit) ...
%                 '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LD_fish' F];
%             dpath = [char(dp) '/'];
%             load([datap dpath sname 'lastyr_sum_mean_biom']);
%             
%             
%             %% Fishing
%             Ftot = sum(Ftotcatch(:,s));
%             Ptot = sum(Ptotcatch(:,s));
%             Dtot = sum(Dtotcatch(:,s));
%             Tot = Ftot+Ptot+Dtot;
%             
%             f7 = figure(7);
%             plot(i,Dtot,'.k','MarkerSize',20); hold on;
%             xlim([0 length(frate)+1])
%             if (i==length(frate))
%                 set(gca,'XTick',1:length(frate),'XTickLabel',fsim);
%                 ylabel('All D')
%                 title(['RE=' num2str(reff(r)) ' ' loc])
%                 stamp(cfile2)
%             end
%         end
%         print(f7,'-dpng',[figp loc '/' sname cfile2 '_' lname 'D.png'])
%         
%     end
    
    %%
    for s=1:length(spots)
        
        loc = spots{s};
        lname = [loc '_'];
        
        %%
        for i=1:length(frate)
            F = frate{i};
            dp = ['Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LD_fish' F];
            dpath = [char(dp) '/'];
            load([datap dpath sname 'lastyr_sum_mean_biom']);
            
            
            %% Fishing
            Ftot = sum(Ftotcatch(:,s));
            Ptot = sum(Ptotcatch(:,s));
            Dtot = sum(Dtotcatch(:,s));
            Tot = Ftot+Ptot+Dtot;
            
            f12 = figure(12);
            subplot(4,3,s)
            plot(fnum(i),Dtot,'.k','MarkerSize',25); hold on;
            xlim([0 fnum(end)+0.1])
            if (i==length(frate))
                %set(gca,'XTick',1:length(frate),'XTickLabel',fsim);
                if (s==2)
                    str = {['RE=' num2str(reff(r))], loc};
                    title(str)
                else
                    title(loc)
                end
                if (s==4)
                    ylabel('Total D catch (g) in final year')
                end
                stamp(cfile2)
            end
            
        end
        
    end
    print(f12,'-dpng',[figp sname cfile2 '_allD_locs.png'])
end
