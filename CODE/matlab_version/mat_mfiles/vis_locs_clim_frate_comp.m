% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Bio_rates/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
sname = 'Climatol_';
H = {'All_fish03','fish_F030_P015_D015','fish_F030_P060_D060',...
    'fish_F015_P030_D030','fish_F060_P030_D030'};
dpath = [datap char(dp) '/'];
cfile = char(dp);
    
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

A = 4.39;
fc = 0.2;
f0 = 0.6;
epsassim = 0.7;
n = 3/4;

AvailEnergy = A*mass.^n;
Consumption = A / (epsassim*(f0-fc)) * mass.^n;

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%%
np = 5;
fishsp = NaN*ones(4,length(spots),np);
mconW = NaN*ones(3,np);
mconF = NaN*ones(length(spots),3,np);
mprod = NaN*ones(3,length(spots),np);
mpred = NaN*ones(3,length(spots),np);
ggeA = NaN*ones(3,length(spots),np);
for j = 1:5
    harv = H{j};
    load([dpath sname 'locs_' harv '_lastyr_sum_mean_biom.mat'])

    %% Biomass of each type
    fishsp(:,:,j) = squeeze(nansum(all_mean));
    
    %% Consump vs. weight
    Pcon = Pcon .* mass .* 365;
    Fcon = Fcon .* mass(1:2,:) .* 365;
    Dcon = Dcon .* mass .* 365;
    Fcon(3,:) = nan;
    Tcon = [Fcon,Pcon,Dcon];
    Mcon = nanmean(Tcon,2);
    mconW(:,j) = Mcon;
    
    %% Consump per biomass (I) by type
    % conF(:,1)=squeeze(nanmean(SF(lyr,8,:)))+squeeze(nanmean(MF(lyr,8,:)));
    % conF(:,2)=squeeze(nanmean(SP(lyr,8,:)))+squeeze(nanmean(MP(lyr,8,:)))+squeeze(nanmean(LP(lyr,8,:)));
    % conF(:,3)=squeeze(nanmean(SD
    mconF(:,:,j) = conF;

    %% Production (= nu * biom)
    % Pprod=[SP_prod;MP_prod;LP_prod];
    % Fprod=[SF_prod;MF_prod];
    % Dprod
    mprod(1,:,j) = Fprod(2,:);
    mprod(2,:,j) = Pprod(3,:);
    mprod(3,:,j) = Dprod(3,:);
    
    %% Predation
    % Ppred=[SP_pred;MP_pred;LP_pred];
    % Fpred=[SF_pred;MF_pred];
    % Dpred=[SD_pred;MD_pred;LD_pred];
    mpred(1,:,j) = nansum(Fpred);
    mpred(2,:,j) = nansum(Ppred);
    mpred(3,:,j) = nansum(Dpred);
    
    %% Gross growth efficiency (= nu/consump)
    % Pgge=[SP_gge;MP_gge;LP_gge];
    % Fgge=[SF_gge;MF_gge];
    % Dgge=[SD_gge;MD_gge;LD_gge];
    ggeA(1,:,j) = Fgge(2,:);
    ggeA(2,:,j) = Pgge(3,:);
    ggeA(3,:,j) = Dgge(3,:);
    
end 
% Save values for all locs 
save([dpath 'Locs_harv_comp.mat'],'fishsp','mconW','mconF','mprod',...
    'mpred','ggeA');

%%
FPrat = squeeze(fishsp(1,:,:)./(fishsp(1,:,:)+fishsp(2,:,:)));
DPrat = squeeze(fishsp(3,:,:)./(fishsp(3,:,:)+fishsp(2,:,:)));

y0 = zeros(7,1);

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f1=figure(1);
    subplot(4,4,s)
    plot(0:6,y0,'--k'); hold on;
    plot(1-0.25:np,log10(squeeze(fishsp(1,s,:)-fishsp(1,s,1))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,log10(squeeze(fishsp(2,s,:)-fishsp(2,s,1))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,log10(squeeze(fishsp(3,s,:)-fishsp(3,s,1))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([-5 1])
    if (s==5)
        ylabel('change in log10 Mean Biom (g m^-^2)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    if (s>12)
        for t=1:np
            text(t,-5.5,H{t},'Rotation',45,'HorizontalAlignment','right')
        end
    end
    stamp(cfile)
    title(loc)
    
    f2=figure(2);
    subplot(4,4,s)
    plot(0:6,y0,'--k'); hold on;
    bar(FPrat(s,:)-FPrat(s,1)); hold on;
    ylim([-0.3 0.3])
    xlim([0 6])
    if (s==2)
        title({'change in Frac F:P in final year'; loc})
    else
        title(loc)
    end
    if (s>12)
        for t=1:np
            text(t,-0.4,H{t},'Rotation',45,'HorizontalAlignment','right')
        end
    end
    stamp(cfile)
    
    %confF
    f3=figure(3);
    subplot(4,4,s)
    plot(0:6,y0,'--k'); hold on;
    plot(1-0.25:np,(squeeze(mconF(s,1,:)-mconF(s,1,1))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mconF(s,2,:)-mconF(s,2,1))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mconF(s,3,:)-mconF(s,3,1))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([-5e-3 2e-3])
    if (s==5)
        ylabel('change in Biom F consumed (g g^-^1 d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    if (s>12)
        for t=1:np
            text(t,-5.2e-3,H{t},'Rotation',45,'HorizontalAlignment','right')
        end
    end
    stamp(cfile)
    title(loc)
    
    %prod
    f4=figure(4);
    subplot(4,4,s)
    plot(0:6,y0,'--k'); hold on;
    plot(1-0.25:np,(squeeze(mprod(1,s,:)-mprod(1,s,1))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mprod(2,s,:)-mprod(2,s,1))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mprod(3,s,:)-mprod(3,s,1))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    %ylim([-2 2])
    if (s==5)
        ylabel('change in adult production (g g^-^1 d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    if (s>12)
        for t=1:np
            text(t,-0.06,H{t},'Rotation',45,'HorizontalAlignment','right')
        end
    end
    stamp(cfile)
    title(loc)
    
    %pred
    f5=figure(5);
    subplot(4,4,s)
    plot(0:6,y0,'--k'); hold on;
    plot(1-0.25:np,(squeeze(mpred(1,s,:)-mpred(1,s,1))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(mpred(2,s,:)-mpred(2,s,1))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(mpred(3,s,:)-mpred(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    %ylim([-2 2])
    if (s==5)
        ylabel('change in predation rate (d^-^1)')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    if (s>12)
        for t=1:np
            text(t,-0.015,H{t},'Rotation',45,'HorizontalAlignment','right')
        end
    end
    stamp(cfile)
    title(loc)
    
    %gge
    f6=figure(6);
    subplot(4,4,s)
    plot(0:6,y0,'--k'); hold on;
    plot(1-0.25:np,(squeeze(ggeA(1,s,:)-ggeA(1,s,1))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,(squeeze(ggeA(2,s,:)-ggeA(2,s,1))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,(squeeze(ggeA(3,s,:)-ggeA(3,s,1))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([-0.2 0.2])
    if (s==5)
        ylabel('change in adult gross growth efficiency')
    end
    set(gca,'XTick',1:np,'XTickLabel',[]);
    if (s>12)
        for t=1:np
            text(t,-0.22,H{t},'Rotation',45,'HorizontalAlignment','right')
        end
    end
    stamp(cfile)
    title(loc)
    
end %spots
print(f1,'-dpng',[figp sname '_harv_comp_tot_mean_biomass_type_all_locs.png'])
print(f2,'-dpng',[figp sname '_harv_comp_FP_frac_all_locs.png'])
print(f3,'-dpng',[figp sname '_harv_comp_conF_all_locs.png'])
print(f4,'-dpng',[figp sname '_harv_comp_Aprod_all_locs.png'])
print(f5,'-dpng',[figp sname '_harv_comp_tpred_all_locs.png'])
print(f6,'-dpng',[figp sname '_harv_comp_Agge_all_locs.png'])

%% Consump vs. weight
figure(7)
subplot(3,1,1)
bar(mconW(1,:))
set(gca,'XTick',1:np,'XTickLabel',[]); hold on;
plot(0:np+1, Consumption(1,1:np+2),'--k','LineWidth',2)
ylabel('S')
title('Mean consumption (g yr^-^1)')

subplot(3,1,2)
bar(mconW(2,:))
set(gca,'XTick',1:np,'XTickLabel',[]); hold on;
plot(0:np+1, Consumption(2,1:np+2),'--k','LineWidth',2)
ylabel('M')

subplot(3,1,3)
bar(mconW(3,:))
set(gca,'XTick',1:np,'XTickLabel',[]); hold on;
for t=1:np
    text(t,-0.1,H{t},'Rotation',45,'HorizontalAlignment','right')
end
plot(0:np+1, Consumption(3,1:np+2),'--k','LineWidth',2)
ylabel('L')
xlabel('Harv')
print('-dpng',[figp sname '_harv_comp_mean_con_size_all_locs.png'])

