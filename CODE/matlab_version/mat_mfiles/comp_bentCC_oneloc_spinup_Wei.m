% Calculate different skill metrics for each oneloc bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
datap = '/Volumes/GFDL/CSV/Matlab_big_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Comparisons/';

RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
CarCap = {'050','100','150','200'};
car = 0.5:0.5:2.0;
benteff = {'05','10'};
beff = 0.05:0.05:0.1;
fcrit = 40;
nmort = '0';
kad = 100;
pref = 'D100';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
sname = 'Spinup_';

seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei_inverts_fish_gWWm2_locs.csv',1,1);

%% colors
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


%% Wei benthic biomass
% in g WW/m2
invert = seafl(:,1);
fish = seafl(:,2);

%% Skill metrics

% Stowe et al 2009
% Univariate (model vs observ)
% 1. correlation coefficient [-1 1]; want close to 1
% 2. root mean square error; want close to 0; considering the magnitude rather than the direction
% 3. average error (bias); want close to 0; neg & pos discrepancies can cancel each other
% 4. average absolute error; want close to 0; considering the magnitude rather than the direction
% 5. modeling efficiency; want close to 1; <0 means avg obs better than model; how well a model predicts relative to avg obs

metrics={'r','RMSE','AE','AAE','MEF'};
metrics=metrics';

%Sometimes it is appropriate to log-transform the observations and predictions before
%calculating goodness-of-fit statistics so that differences between predicted and observed
%values will not be highly skewed and dominated by a small proportion of high values.

ndp = length(CarCap);

r = NaN*ones(length(RE),ndp,length(benteff));
RMSE = r;
AE = r;
AAE = r;
MEF = r;

for i=1:length(benteff)
    BE = benteff{i};
    for R = 1:length(RE)
        rfrac = RE{R};
        
        cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_RE' rfrac...
                '_BentCCtests'];
        load([datap 'Bent_CC_tests/' cfile2 '.mat'],'Bbio','MLD');
        
        close all
        
        %% Model bent & dem
        modi = Bbio';
        modf = MLD';
        
        %% Bar plots for comparison
        for s=1:length(spots)
            loc = spots{s};
            lname = [loc '_'];
            
            f1=figure(1);
            subplot(4,3,s)
            bar(modi(s,:)); hold on;
            plot(0:ndp+1,repmat(invert(s),ndp+2,1),'r','LineWidth',2)
            xlim([0 ndp+1])
            if (s==4)
                ylabel('Mean Biom (g m^-^2) in final year')
            end
            set(gca,'XTick',1:ndp,'XTickLabel',car);
            if (s==1 || s==2 || s==6 || s==10)
                ylim([0 0.6])
            else
                ylim([0 0.1])
            end
            if (s==11)
                text(8,0.1,['BE=' num2str(beff(i))]);
                text(8,0,['RE=' num2str(reff(R))]);
            end
            stamp(cfile2)
            title([loc ' B biom'])
            
            f2=figure(2);
            subplot(4,3,s)
            bar(log10(modf(s,:))); hold on;
            plot(0:ndp+1,repmat(log10(fish(s)),ndp+2,1),'r','LineWidth',2)
            xlim([0 ndp+1])
            ylim([-3 1])
            if (s==4)
                ylabel('log10 Mean Biom (g m^-^2) in final year')
            end
            set(gca,'XTick',1:ndp,'XTickLabel',car);
            if (s==11)
                text(8,0,['BE=' num2str(beff(i))]);
                text(8,-1,['RE=' num2str(reff(R))]);
            end
            stamp(cfile2)
            title([loc ' D biom'])
        end
        print(f1,'-dpng',[figp sname cfile2 '_Wei_bent_biomass_comp.png'])
        %print(f2,'-dpng',[figp sname cfile2 '_Wei_fish_biomass_comp.png'])
        
        
        %% Inverts
        obs = real(log(invert));
        model = real(log(modi));
        
        n = length(obs);
        
        iskill=NaN*ones(5,ndp);
        for j=1:ndp
            o=obs;
            p=model(:,j);
            p(isnan(p))=0;
            
            omean=repmat(nanmean(o),n,1);
            pmean=repmat(nanmean(p),n,1);
            osig=nanstd(o);
            psig=nanstd(p);
            
            % corr coeff
            num=nansum((o-omean).*(p-pmean));
            d1=nansum((o-omean).^2);
            d2=nansum((p-pmean).^2);
            den=sqrt(d1*d2);
            iskill(1,j) = num/den;
            
            % root mean square error
            num=nansum((p-o).^2);
            iskill(2,j) = sqrt(num/n);
            
            % average error
            iskill(3,j) = nansum(p-o) / n;
            
            % average absolute error
            iskill(4,j) = nansum(abs(p-o)) / n;
            
            % modeling efficiency
            num1=nansum((o-omean).^2);
            num2=nansum((p-o).^2);
            iskill(5,j) = (num1-num2)/num1;
        end
        
        %% Fish
        obs = real(log(fish));
        model = real(log(modf));
        
        n = length(obs);
        
        fskill=NaN*ones(5,ndp);
        for j=1:ndp
            o=obs;
            p=model(:,j);
            p(isnan(p))=0;
            p(isinf(p))=0;
            
            omean=repmat(nanmean(o),n,1);
            pmean=repmat(nanmean(p),n,1);
            osig=nanstd(o);
            psig=nanstd(p);
            
            % corr coeff
            num=nansum((o-omean).*(p-pmean));
            d1=nansum((o-omean).^2);
            d2=nansum((p-pmean).^2);
            den=sqrt(d1*d2);
            fskill(1,j) = num/den;
            
            % root mean square error
            num=nansum((p-o).^2);
            fskill(2,j) = sqrt(num/n);
            
            % average error
            fskill(3,j) = nansum(p-o) / n;
            
            % average absolute error
            fskill(4,j) = nansum(abs(p-o)) / n;
            
            % modeling efficiency
            num1=nansum((o-omean).^2);
            num2=nansum((p-o).^2);
            fskill(5,j) = (num1-num2)/num1;
        end %CC
        
        %% Plot results
        
        % Bar graphs
        figure(4)
        subplot(2,2,1)
        bar(iskill(1,:),'k')
        ylabel('Correlation coefficient')
        ylim([-0.1 1])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
        stamp(cfile2)
        title(['Inverts skill with BE=' num2str(beff(i)) ' and RE=' num2str(reff(R))]...
            ,'HorizontalAlignment','left')
        
        subplot(2,2,2)
        bar(iskill(2,:),'k')
        ylabel('Root mean square error')
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        ylim([0 2.5])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
        
        subplot(2,2,3)
        bar(iskill(3,:),'k')
        ylabel('Average error')
        ylim([0 2.5])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
        
        subplot(2,2,4)
        bar(iskill(5,:),'k')
        ylabel('Modeling efficiency')
        ylim([-5 1])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
        print('-dpng',[figp sname cfile2 '_locs_skill_Wei_inverts'])
        
        %%
        figure(5)
        subplot(2,2,1)
        bar(fskill(1,:),'k')
        ylabel('Correlation coefficient')
        ylim([-0.1 1])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
        stamp(cfile2)
        title(['Fish skill with BE=' num2str(beff(i)) ' and RE=' num2str(reff(R))]...
            ,'HorizontalAlignment','left')
        
        subplot(2,2,2)
        bar(fskill(2,:),'k')
        ylabel('Root mean square error')
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        %ylim([0 2.5])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-0.1,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
%         
        subplot(2,2,3)
        bar(fskill(3,:),'k')
        ylabel('Average error')
        %ylim([-1 1])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-1.1,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
%         
        subplot(2,2,4)
        bar(fskill(5,:),'k')
        ylabel('Modeling efficiency')
        %ylim([-1 1])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
%         for t=2:2:ndp
%             text(t,-1.1,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
%         end
        %print('-dpng',[figp sname cfile2 '_locs_skill_Wei_fish'])
        
        %% Best
        one=[1;5];
        zer=[2;3;4];
        for z=1:5
            Z1=sum(z==one);
            Z2=sum(z==zer);
            if (Z1>0)
                ibest{z,1}=car(find(min(abs(1-(iskill(z,:))))==abs(1-(iskill(z,:)))));
                fbest{z,1}=car(find(min(abs(1-(fskill(z,:))))==abs(1-(fskill(z,:)))));
            elseif (Z2>0)
                ibest{z,1}=car(find(min(abs(0-(iskill(z,:))))==abs(0-(iskill(z,:)))));
                fbest{z,1}=car(find(min(abs(0-(fskill(z,:))))==abs(0-(fskill(z,:)))));
            end
        end
        
        Ti=table(metrics,ibest);
        Tf=table(metrics,fbest);
        
        save([datap 'Bent_CC_tests/' cfile2 '_locs_skill_Wei.mat'],'modi',...
            'modf','invert','fish','metrics','iskill','fskill','Ti','Tf');
        
        r(R,:,i) = iskill(1,:);
        RMSE(R,:,i) = iskill(2,:);
        AE(R,:,i) = iskill(3,:);
        AAE(R,:,i) = iskill(4,:);
        MEF(R,:,i) = iskill(5,:);
        
    end %RE
end %BE

%% Find best
[Cmesh,Rmesh,Bmesh] = meshgrid(car,reff,beff);

maxr = find(r==max(r(:)));
maxMEF = find(MEF==max(MEF(:)));
zerAE = find( abs(0-(AE(:))) == min(abs(0-(AE(:)))) );
zerAAE = find( abs(0-(AAE(:))) == min(abs(0-(AAE(:)))) );
zerRMSE = find( abs(0-(RMSE(:))) == min(abs(0-(RMSE(:)))) );

value = [max(r(:)); min(abs(0-(RMSE(:)))); min(abs(0-(AE(:)))); min(abs(0-(AAE(:))));...
    max(MEF(:))];

Ibest(1,1) = Bmesh(maxr);
Ibest(1,2) = Rmesh(maxr);
Ibest(1,3) = Cmesh(maxr);

Ibest(2,1) = Bmesh(zerRMSE);
Ibest(2,2) = Rmesh(zerRMSE);
Ibest(2,3) = Cmesh(zerRMSE);

Ibest(3,1) = Bmesh(zerAE);
Ibest(3,2) = Rmesh(zerAE);
Ibest(3,3) = Cmesh(zerAE);

Ibest(4,1) = Bmesh(zerAAE);
Ibest(4,2) = Rmesh(zerAAE);
Ibest(4,3) = Cmesh(zerAAE);

Ibest(5,1) = Bmesh(maxMEF);
Ibest(5,2) = Rmesh(maxMEF);
Ibest(5,3) = Cmesh(maxMEF);

I=table(metrics,value,Ibest(:,1),Ibest(:,2),Ibest(:,3),'VariableNames',...
    {'Metric','Value','BentEff','RepEff','CarCap'});

cfile3 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort'  nmort '_BentCCtests'];
save([datap 'Bent_CC_tests/' cfile3 '_locs_iskill_Wei.mat'],'r',...
            'RMSE','AE','AAE','MEF','beff','reff','car','I');
writetable(I,[datap 'Bent_CC_tests/' cfile3 '_locs_best_iskill_Wei.csv'])
        
