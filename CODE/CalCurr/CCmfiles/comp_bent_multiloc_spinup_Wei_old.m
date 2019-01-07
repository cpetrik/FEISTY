% Calculate different skill metrics for each oneloc bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Mat_runs/';
datap = '/Volumes/GFDL/CSV/Matlab_orig_size/';

RE = {'1000','0500','0100'};%,'0050','0010','00050','00010'};
reff = [1.0,0.5,0.1];%,0.05,0.01,0.005,0.001];
CarCap = {'050','100','150','200','250','300'};
car = 0.5:0.5:3.0;
benteff = {'05','10'};%,'15','20'};
beff = 0.05:0.05:0.1;
fcrit = 40;
nmort = '0';
kad = 100;
pref = 'D100';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};

sname = 'Spinup_';

%% Wei benthic biomass
seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei2010_Global_seafloor_biomass.csv',1,0);
Wcol = {'latitude','longitude','depth','bact.biom.mean','meio.biom.mean',...
    'macro.biom.mean','mega.biom.mean','inv.biom.mean','fis.biom.mean'};
Wcol = Wcol';

% all mean biomasses in log10 mg C/m2
invert = seafl(:,8);
fish = seafl(:,9);
% convert to g WW/m2
invert = 10.^(invert) * 1e-3 * 9.0;
fish = 10.^(fish) * 1e-3 * 9.0;

% put on same grid as POEM output
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

Zi=griddata(seafl(:,2),seafl(:,1),invert,X,Y);
Zf=griddata(seafl(:,2),seafl(:,1),fish,X,Y);

%% Maps
% Inv
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zi))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Wei benthic invert biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_inverts.png')

% Fish
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Wei benthic fish biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_fish.png')

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
NSTD = r;
RMSD = r;

%%
for i=1:length(benteff)
    BE = benteff{i};
    for R = 1:length(RE)
        rfrac = RE{R};
        iskill=NaN*ones(7,ndp);
        fskill=NaN*ones(7,ndp);
        ibest=cell(5,1);
        fbest=cell(5,1);
        for c = 1:length(CarCap)
            CC = CarCap{c};
            
            %%
            cfile = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_CC' CC '_RE' rfrac];
            fpath=['/Volumes/GFDL/NC/Matlab_runs/' cfile '/'];
            ppath = [pp cfile];
            load([fpath 'Means_spinup_' cfile '.mat'],'md_mean','ld_mean','b_mean');
            
            cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_RE' rfrac '_BentCCtests'];
            
            close all
            
            %% Model bent & dem
            Zmd=griddata(grid(:,2),grid(:,3),md_mean,X,Y);
            Zld=griddata(grid(:,2),grid(:,3),ld_mean,X,Y);
            Zb=griddata(grid(:,2),grid(:,3),b_mean,X,Y);
            
            Zd = Zmd+Zld;
            
            % Bent
            figure(1)
            m_proj('miller','lat',82);
            m_pcolor(X,Y,real(log10(Zb))); hold on;
            shading flat
            m_coast('patch',[.5 .5 .5],'edgecolor','none');
            m_grid;
            title('log10 mean model benthic biomass (g m^-^2)')
            colormap('jet')
            colorbar('h')
            caxis([-3 1])
            stamp(cfile)
            print('-dpng',[ppath '_Spinup_global_Bent.png'])
            
            % D
            figure(2)
            m_proj('miller','lat',82);
            m_pcolor(X,Y,real(log10(Zd))); hold on;
            shading flat
            m_coast('patch',[.5 .5 .5],'edgecolor','none');
            m_grid;
            title('log10 mean model M&L D biomass (g m^-^2)')
            colormap('jet')
            colorbar('h')
            caxis([-3 1])
            stamp(cfile)
            print('-dpng',[ppath '_Spinup_global_Demersal.png'])
            
            
            %% Inverts
            obs = real(log(Zi(:)));
            model = real(log(Zb(:)));
            
            n = length(obs);
            
            o=obs;
            p=model;
            %             o(isnan(o))=0;
            %             p(isnan(p))=0;
            
            omean=repmat(nanmean(o),n,1);
            pmean=repmat(nanmean(p),n,1);
            osig=nanstd(o);
            psig=nanstd(p);
            
            % corr coeff
            num=nansum((o-omean).*(p-pmean));
            d1=nansum((o-omean).^2);
            d2=nansum((p-pmean).^2);
            den=sqrt(d1*d2);
            iskill(1,c) = num/den;
            
            % root mean square error
            num=nansum((p-o).^2);
            iskill(2,c) = sqrt(num/n);
            
            % average error
            iskill(3,c) = nansum(p-o) / n;
            
            % average absolute error
            iskill(4,c) = nansum(abs(p-o)) / n;
            
            % modeling efficiency
            num1=nansum((o-omean).^2);
            num2=nansum((p-o).^2);
            iskill(5,c) = (num1-num2)/num1;
            
            % Taylor normalized std
            iskill(6,c) = psig/osig;
            
            % total root mean square difference
            q1=nansum((o-p).^2);
            iskill(7,c) = sqrt(q1/n);
            
            
            %% Fish
            obs = real(log(Zf(:)));
            model = real(log(Zd(:)));
            
            n = length(obs);
            
            o=obs;
            p=model;
            
            omean=repmat(nanmean(o),n,1);
            pmean=repmat(nanmean(p),n,1);
            osig=nanstd(o);
            psig=nanstd(p);
            
            % corr coeff
            num=nansum((o-omean).*(p-pmean));
            d1=nansum((o-omean).^2);
            d2=nansum((p-pmean).^2);
            den=sqrt(d1*d2);
            fskill(1,c) = num/den;
            
            % root mean square error
            num=nansum((p-o).^2);
            fskill(2,c) = sqrt(num/n);
            
            % average error
            fskill(3,c) = nansum(p-o) / n;
            
            % average absolute error
            fskill(4,c) = nansum(abs(p-o)) / n;
            
            % modeling efficiency
            num1=nansum((o-omean).^2);
            num2=nansum((p-o).^2);
            fskill(5,c) = (num1-num2)/num1;
            
            % Taylor normalized std
            fskill(6,c) = psig/osig;
            
            % total root mean square difference
            q1=nansum((o-p).^2);
            fskill(7,c) = sqrt(q1/n);
            
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
        ylim([0 1])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
        %         for t=2:2:ndp
        %             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
        %         end
        
        subplot(2,2,3)
        bar(iskill(3,:),'k')
        ylabel('Average error')
        ylim([-0.1 0.4])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
        %         for t=2:2:ndp
        %             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
        %         end
        
        subplot(2,2,4)
        bar(iskill(5,:),'k')
        ylabel('Modeling efficiency')
        ylim([0 0.5])
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',car);
        %         for t=2:2:ndp
        %             text(t,-0.15,num2str(car(t)),'Rotation',45,'HorizontalAlignment','right')
        %         end
        print('-dpng',[ppath '_Spinup_skill_Wei_inverts'])
        
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
        print('-dpng',[ppath '_Spinup_skill_Wei_fish'])
        
        %% Taylor diagrams
        % INVERTS
        [rmsd,it]=sort(iskill(7,:),'descend');
        theta=acos(iskill(1,it));    %corr coeff
        rho=iskill(6,it);            %stdev
        simtex=CarCap(it);
        simtex{ndp+1}='obs';
        
        % % Just those <2
        % id = find(rho<2);
        % theta = theta(id);
        % rho = rho(id);
        % simtex=sims(id);
        % simtex{length(id)+1}='obs';
        
        tr=0;
        rr=1;
        figure(6)
        h0=polar(0,1.5,'.'); hold on;
        set(h0,'color','w');
        for s=1:ndp
            h=polar(theta(s),rho(s),'.'); hold on;
            set(h,'color',cm{s},'MarkerSize',25);
            %     set(h,'MarkerSize',25);
        end
        h2=polar(tr,rr,'k*');
        set(h2,'MarkerSize',10);
        axis([0 1.5 0 1.5])
        title('Inverts Taylor diagram')
        legend([' ' simtex])
        legend('location','northeast')
        print('-dpng',[ppath '_Spinup_global_invert_Taylor_Wei'])
        
        %% FISH
        
        [rmsd,it]=sort(fskill(7,:),'descend');
        theta=acos(fskill(1,it));    %corr coeff
        rho=fskill(6,it);            %stdev
        simtex=CarCap(it);
        simtex{ndp+1}='obs';
        
        % Just those <1.5
        % id = find(rho<1.5);
        % theta = theta(id);
        % rho = rho(id);
        % simtex=sims(id);
        % simtex{length(id)+1}='obs';
        
        tr=0;
        rr=1;
        figure(7)
        h0=polar(0,1.5,'.'); hold on;
        set(h0,'color','w');
        for s=1:ndp
            h=polar(theta(s),rho(s),'.'); hold on;
            set(h,'color',cm{s},'MarkerSize',25);
            %     set(h,'MarkerSize',25);
        end
        h2=polar(tr,rr,'k*');
        set(h2,'MarkerSize',10);
        axis([0 1.5 0 1.5])
        title('Fish Taylor diagram')
        legend([' ' simtex])
        legend('location','northeast')
        print('-dpng',[ppath '_Spinup_global_fish_Taylor_Wei'])
        
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
        
        save([datap 'Bent_CC_tests/' cfile2 '_global_skill_Wei.mat'],...
            'invert','fish','metrics','iskill','fskill','Ti','Tf');
        
        r(R,:,i) = iskill(1,:);
        RMSE(R,:,i) = iskill(2,:);
        AE(R,:,i) = iskill(3,:);
        AAE(R,:,i) = iskill(4,:);
        MEF(R,:,i) = iskill(5,:);
        NSTD(R,:,i) = iskill(6,:);
        RMSD(R,:,i) = iskill(7,:);
        
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

%%
cfile3 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
    '_' pref '_nmort'  nmort '_BentCCtests'];
save([datap 'Bent_CC_tests/' cfile3 '_global_iskill_Wei.mat'],'r',...
    'RMSE','AE','AAE','MEF','beff','reff','car','I');
writetable(I,[datap 'Bent_CC_tests/' cfile3 '_global_best_iskill_Wei.csv'])

