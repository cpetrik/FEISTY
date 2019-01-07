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
encs = linspace(10,100,10); %logspace(2,4,10); %
cmaxs = linspace(10,100,10);
efn = 70;
tefn = num2str(efn);
cfn = 50;
tcfn = num2str(cfn);
fcrit = 20;
nmort = '1';
kad = 50;
D = 'D075';
J = 'J075';
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

A = 4.39;
fc = 0.2;
f0 = 0.6;
epsassim = 0.7;
n = 3/4;

AvailEnergy = A*mass.^n;
Consumption = A / (epsassim*(f0-fc)) * mass.^n;

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%%
for R = 3;%1:length(RE)
    rfrac = RE{R};
    close all
    
        fishsp  = NaN*ones(4,length(spots),length(encs));
        mcon  = NaN*ones(3,length(encs));
        mlev  = NaN*ones(3,length(encs));
        for E = 1:length(encs)
            efn = encs(E);
            tefn = num2str(round(efn));
            dp = ['Dc_enc',tefn,'_cmax-metab',tcfn,'_fcrit',num2str(fcrit),'_',D,'_',J,'_',Sm,...
                '_nmort',nmort,'_BE',BE,'_CC',CC,'_RE',rfrac];
            dpath = [datap char(dp) '/'];
            fpath = [figp char(dp) '/'];
            cfile = char(dp);
            cfile2 = ['Dc_enc',tefn,'_cmax-metab_fcrit',num2str(fcrit),'_',D,'_',J,'_',Sm,...
                '_nmort',nmort,'_BE',BE,'_CC',CC,'_RE',rfrac,'_Enctests'];
    
            load([dpath sname 'lastyr_sum_mean_biom.mat']);
    
            %% Biomass of each type
            fishsp(:,:,E) = squeeze(nansum(all_mean));
    
            %% Consump vs. weight
            Pcon = Pcon .* mass .* 365;
            Fcon = Fcon .* mass(1:2,:) .* 365;
            Dcon = Dcon .* mass .* 365;
            Fcon(3,:) = nan;
            Tcon = [Fcon,Pcon,Dcon];
            Mcon = nanmean(Tcon,2);
            mcon(:,E) = Mcon;
    
            %% Feeding level
            Flev(3,:) = nan;
            Tlev = [Flev,Plev,Dlev];
            Mlev = nanmean(Tlev,2);
            mlev(:,E) = Mlev;
    
        end %E
end %RE

%% Save values for all locs and all RE that combo
save([datap 'Bio_rates/' cfile2 '.mat'],'fishsp','mcon','mlev');

%%
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f2=figure(2);
    subplot(4,3,s)
    plot(1-0.25:10,log10(squeeze(fishsp(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:10,log10(squeeze(fishsp(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:11,log10(squeeze(fishsp(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 11])
    ylim([-2 2])
    if (s==4)
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    set(gca,'XTick',1:10,'XTickLabel',[]);
    for t=1:10
        %text(t,-2.1,num2str(cmaxs(t)),'Rotation',45,'HorizontalAlignment','right')
        text(t,-2.1,num2str(encs(t)),'Rotation',45,'HorizontalAlignment','right')
    end
    if (s==11)
        text(15,1,['D=' D]);
        text(15,0.5,['J=' J]);
        text(15,0,['Sm=' Sm]);
        %text(15,-0.5,['enc=' tefn]);
        text(15,-0.5,['cmax-metab=' tcfn]);
    end
    stamp(cfile2)
    title([loc ' All stages'])
    
end %spots
print(f2,'-dpng',[figp sname cfile2 '_tot_mean_biomass_type_all_locs.png'])

%% Encs
% Consump vs. weight
figure(1)
subplot(3,1,1)
bar(mcon(1,:))
set(gca,'XTick',1:10,'XTickLabel',encs); hold on;
plot(0.5:10.5, Consumption(1,:),'--k','LineWidth',2)
ylabel('S')
title('Mean consumption (g yr^-^1)')

subplot(3,1,2)
bar(mcon(2,:))
set(gca,'XTick',1:10,'XTickLabel',encs); hold on;
plot(0.5:10.5, Consumption(2,:),'--k','LineWidth',2)
ylabel('M')

subplot(3,1,3)
bar(mcon(3,:))
set(gca,'XTick',1:10,'XTickLabel',encs); hold on;
plot(0.5:10.5, Consumption(3,:),'--k','LineWidth',2)
ylabel('L')
xlabel('Enc coeff')
print('-dpng',[figp sname cfile2 '_mean_con_size_all_locs.png'])

% Feeding level
figure(3)
subplot(3,1,1)
bar(mlev(1,:))
ylim([0 1])
set(gca,'XTick',1:10,'XTickLabel',encs,'YTick',0.2:0.2:0.8)
ylabel('S')
title('Mean feeding level')

subplot(3,1,2)
bar(mlev(2,:))
ylim([0 1])
set(gca,'XTick',1:10,'XTickLabel',encs,'YTick',0.2:0.2:0.8); 
ylabel('M')

subplot(3,1,3)
bar(mlev(3,:))
ylim([0 1])
set(gca,'XTick',1:10,'XTickLabel',encs,'YTick',0.2:0.2:0.8); 
ylabel('L')
xlabel('Enc coeff')
print('-dpng',[figp sname cfile2 '_mean_flev_size_all_locs.png'])





