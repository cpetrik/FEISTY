% Visualize output of POEM MF fishing tests
% Spinup at one location
% 50 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Fishing/Ponly_';

fnum = [0.1:0.1:1];%,1.2:0.2:2];
fsim = {'.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};%,'1.2','1.4','1.6','1.8','2'};
frate = {'01','02','03','04','05','06','07','08','09','10'};%,'12','14','16','18','20'};
RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
% nmrt = [0,2:5];
% Mort = {'None','Hartvig','mizer','J&C','P&W'};
nmort = '0';
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
    for i=1:length(frate)
        F = frate{i};
        dp = ['PonlyDc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LP_fish' F];
        dpath = [datap char(dp) '/'];
        %fpath = [figp char(dp) '/'];
        fpath = figp;
        cfile = char(dp);
        cfile2 = ['PonlyDc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac ...
            '_LP_fishing_catch'];
        
        %%
        Psum = NaN*ones(3,length(spots));
        Pmean = Psum;
        Pfish = Psum;
        Ptotcatch = Psum;
        all_mean=NaN*ones(3,4,length(spots));
        
        %%
        for s=11;%1:length(spots)
            %%
            loc = spots{s};
            lname = [loc '_'];
            SP = csvread([dpath sname lname 'Sml_p.csv']);
            MP = csvread([dpath sname lname 'Med_p.csv']);
            LP = csvread([dpath sname lname 'Lrg_p.csv']);
            CO = csvread([dpath sname lname 'Cobalt.csv']);
            
            t=1:length(SP);
            lyr=t((end-12+1):end);
            
            %% Final mean biomass in each size
            
            SP_mean=mean(SP(lyr,1));
            MP_mean=mean(MP(lyr,1));
            LP_mean=mean(LP(lyr,1));
            
            P_mean=[SP_mean;MP_mean;LP_mean];
            
            Pmean(:,s) = P_mean;
            
            all_mean(:,2,s) = P_mean;
            
            %% Fishing
            MP_fish=nanmean(MP(lyr,25));
            LP_fish=nanmean(LP(lyr,25));
            
            P_fish=[0;MP_fish;LP_fish];
            
            Pfish(:,s) = P_fish;
            
            MP_totcatch=nansum(MP(lyr,25));
            LP_totcatch=nansum(LP(lyr,25));
            
            P_totcatch=[0;MP_totcatch;LP_totcatch];
            
            Ptotcatch(:,s) = P_totcatch;
            
        end
        
        save([dpath sname 'lastyr_sum_mean_biom'],'Pmean','all_mean',...
            'Pfish','Ptotcatch');
        
    end
    
%     %%
%     for s=1:length(spots)
%         
%         loc = spots{s};
%         lname = [loc '_'];
%         close all
%         %%
%         for i=1:length(frate)
%             F = frate{i};
%             dp = ['Dc_TrefO_cmax-metab2_enc1_MFeqMP_fcrit' num2str(fcrit) ...
%             '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LP_fish' F];
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
%             plot(i,Ptot,'.k','MarkerSize',20); hold on;
%             xlim([0 length(frate)+1])
%             if (i==length(frate))
%                 set(gca,'XTick',1:length(frate),'XTickLabel',fsim);
%                 ylabel('All F')
%                 title(['RE=' num2str(reff(r)) ' ' loc])
%                 stamp(cfile2)
%             end
%         end
%         print(f7,'-dpng',[figp loc '/' sname cfile2 '_' lname 'P.png'])
%         
%     end
    
    %%
    for s=1:length(spots)
        
        loc = spots{s};
        lname = [loc '_'];
        
        %%
        for i=1:length(frate)
            F = frate{i};
            dp = ['PonlyDc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LP_fish' F];
            dpath = [char(dp) '/'];
            load([datap dpath sname 'lastyr_sum_mean_biom']);
            
            
            %% Fishing
            Ptot = sum(Ptotcatch(:,s));
            
            f12 = figure(12);
            subplot(4,3,s)
            plot(fnum(i),Ptot,'.k','MarkerSize',25); hold on;
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
                    ylabel('Total P catch (g) in final year')
                end
                stamp(cfile2)
            end
            
        end
        
    end
    print(f12,'-dpng',[figp sname cfile2 '_allP_locs.png'])
end
