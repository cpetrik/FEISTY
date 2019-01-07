% Visualize output of POEM mortality tests
% Spinup at one location
% 50 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_big_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Feeding/';

RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
nmrt = [0,2:3]; %[0,2:5];
Mort = {'None','Hartvig','mizer'};%,'J&C','P&W'};
fc = 10:10:50;
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
stages = {'SF','MF','SP','MP','LP','SD','MD','LD'};

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

%%
for f = 1:length(fc)
    fcrit = fc(f);
    for i=1:length(nmrt)
        nmort = num2str(nmrt(i));
        for R = 1:length(RE)
            rfrac = RE{R};
            
            dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_CC' CC '_RE' rfrac];
            dpath = [datap char(dp) '/'];
            fpath = [figp char(dp) '/'];
            cfile = char(dp);
            cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_RE' rfrac...
                '_feeding_level'];
            
            load([dpath sname 'lastyr_sum_mean_biom.mat']);
            Lev = [Flev;Plev;Dlev];
            
            %%
            for s=1:length(spots)
                loc = spots{s};
                lname = [loc '_'];
                
                %% Feeding level by stage
                f2=figure(2);
                subplot(4,3,s)
                bar(Lev(:,s))
                xlim([0.5 8.5])
                ylim([0 1])
                if (s==4)
                    ylabel('Mean feeding level in final year')
                end
                set(gca,'XTick',1:length(stages),'XTickLabel',[]);
                for t=1:length(stages)
                    text(t,-0.05,stages{t},'Rotation',45,'HorizontalAlignment','right')
                end
                if (s==11)
                    text(11,0.7,['fcrit=' num2str(fcrit)]);
                    text(11,0.4,['nmort=' nmort]);
                    text(11,0.1,['RE=' num2str(reff(R))]);
                end
                stamp(cfile2)
                title(loc)
                
            end %spots
            print(f2,'-dpng',[figp sname cfile2 '_mean_feeding_lev_all_locs.png'])
            
            % Save values for all locs and all CC for that RE and BE combo
            %save([datap cfile2 '.mat'],'fishsp')
            
        end %RE
    end %nmort
end %fcrit

