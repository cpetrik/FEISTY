% Visualize output of POEM mortality tests
% Spinup at one location
% 50 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_big_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Fishing/';

Frate = 0.1:0.1:1;
fsim = {'.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
frate = {'01','02','03','04','05','06','07','08','09','10'};
RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
nmrt = 0;
% nmrt = [0,2:5];
% Mort = {'None','Hartvig','mizer','J&C','P&W'};
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
nmort = num2str(nmrt);
ndp = length(frate);

for r = 1:length(RE)
    close all
    rfrac = RE{r};
    fishsp  = NaN*ones(4,length(spots),ndp);
    for i=1:length(frate)
        F = frate{i};
        dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LD_fish' F];
        dpath = [datap char(dp) '/'];
        %fpath = [figp char(dp) '/'];
        fpath = figp;
        cfile = char(dp);
        cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac ...
            '_LD_fishing'];
        
        load([dpath sname 'lastyr_sum_mean_biom.mat']);
        
        %% Sum mean biom over stages
        fishsp(:,:,i) = squeeze(nansum(all_mean));
        
    end %Frate
    
    %%
    for s=1:length(spots)
        loc = spots{s};
        lname = [loc '_'];
        
        %% Biomass of all fish in each type
        f2=figure(2);
        subplot(4,3,s)
        plot(1-0.25:ndp,log10(squeeze(fishsp(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:ndp,log10(squeeze(fishsp(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1+0.25:ndp+0.25,log10(squeeze(fishsp(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-2 1.5])
        if (s==4)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,-2.1,fsim{t},'Rotation',45,'HorizontalAlignment','right')
        end
        if (s==11)
            text(14,1,['nmort=' nmort]);
            text(14,0,'Fishing LD');
            text(14,-1,['RE=' num2str(reff(r))]);
        end
        stamp(cfile2)
        %title([loc ' All stages'])
        title(loc)
        
    end %spots
    print(f2,'-dpng',[figp sname cfile2 '_tot_mean_biomass_type_all_locs.png'])
    
    % Save values for all locs and all CC for that RE and BE combo
    %save([datap cfile2 '.mat'],'fishsp')
    
end %nmort

