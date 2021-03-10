% Visualize output of FEISTY Climatology at single locations
% 150 years, monthly means saved
% Test diff values of consumption half saturation
% to get appropriate turnover rate of small fish

clear all
close all

datap = '/Volumes/MIP/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/Climatol/half_sat_exper/'];

%%
load('cmap_ppt_angles.mat')
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
kays = [1.25:0.25:3 3.5:0.5:7 10:20];
%kays = [2:7 10:20];
np = length(kays);
% need: nu, gam, rep, rec, yield, die, mort, con, bio
fishsp  = NaN*ones(4,length(spots),np);
mnu     = NaN*ones(3,length(spots),np);
mgam    = NaN*ones(3,length(spots),np);
mrec    = NaN*ones(3,length(spots),np);
mrep    = NaN*ones(3,length(spots),np);
myield  = NaN*ones(3,length(spots),np);
mdie    = NaN*ones(3,length(spots),np);
mmort   = NaN*ones(3,length(spots),np);
mcon    = NaN*ones(3,length(spots),np);
mlev    = NaN*ones(3,length(spots),np);
mzfrac  = NaN*ones(3,length(spots),np);

%%
for n = 1:length(kays)
    Ka = kays(n);
    tk = num2str(int64(100*Ka));
    close all
    sname = ['Climatol_ksat',tk,'_locs_'];
    cfile = char(dp);
    
    load([dpath sname harv '_lastyr_sum_mean_biom.mat']);
    
    %% Biomass of each type
    fishsp(:,:,n) = squeeze(nansum(all_mean));
    %% Nu - energy for biomass production
    mnu(1,:,n) = Fmgr(1,:);
    mnu(2,:,n) = Pmgr(1,:);
    mnu(3,:,n) = Dmgr(1,:);
    %% Recruitment into stage (make per biomass)
    mrec(1,:,n) = Frec(1,:)./Fmean(1,:);
    mrec(2,:,n) = Prec(1,:)./Pmean(1,:);
    mrec(3,:,n) = Drec(1,:)./Dmean(1,:);
    %% Reproduction rate (per biomass)
    mrep(1,:,n) = Frep(1,:);
    mrep(2,:,n) = Prep(1,:);
    mrep(3,:,n) = Drep(1,:);
    %% Gamma - growth out of stage (per biomass)
    mgam(1,:,n) = Fgam(1,:);
    mgam(2,:,n) = Pgam(1,:);
    mgam(3,:,n) = Dgam(1,:);
    %% Con (per biomass)
    mcon(1,:,n) = Fcon(1,:);
    mcon(2,:,n) = Pcon(1,:);
    mcon(3,:,n) = Dcon(1,:);
    %% Die - predation (per biomass)
    mdie(1,:,n) = Fpred(1,:);
    mdie(2,:,n) = Ppred(1,:);
    mdie(3,:,n) = Dpred(1,:);
    %% Nmort (per biomass)
    mmort(1,:,n) = Fnat(1,:);
    mmort(2,:,n) = Pnat(1,:);
    mmort(3,:,n) = Dnat(1,:);
    %% Fishing rate
    myield(1,:,n) = Ffrate(2,:);
    myield(2,:,n) = Pfrate(3,:);
    myield(3,:,n) = Dfrate(3,:);
    %% Feeding level
    mlev(1,:,n) = Flev(1,:);
    mlev(2,:,n) = Plev(1,:);
    mlev(3,:,n) = Dlev(1,:);
    %% Zoop frac con
    mzfrac(:,:,n) = z';
    
end
% Save values for all locs and all ksats that combo
save([dpath 'Climatol_locs_ksat_125_2000_small.mat'],'fishsp','mnu','mrec',...
    'mrep','mcon','mdie','mmort','myield','mlev','mgam');


%% Figures
mzfrac = min(mzfrac,1);
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Summed mean biom over stages
    f1=figure(1);
    subplot(4,4,s)
    plot(kays,log10(squeeze(fishsp(1,s,:))),'color',cmap_ppt(3,:),'LineWidth',1.5); hold on;
    plot(kays,log10(squeeze(fishsp(2,s,:))),'color',cmap_ppt(1,:),'LineWidth',1.5); hold on;
    plot(kays,log10(squeeze(fishsp(3,s,:))),'color',cmap_ppt(2,:),'LineWidth',1.5); hold on;
    xlim([kays(1) kays(end)])
    ylim([-2 2])
    if (s==5)
        ylabel('log10 Mean Biom (g m^-^2)')
    end
    stamp('')
    title([loc ' All stages'])
    
    %% Nu
    f2=figure(2);
    subplot(4,4,s)
    plot(kays,(squeeze(mnu(1,s,:))),'color',cmap_ppt(3,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mnu(2,s,:))),'color',cmap_ppt(1,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mnu(3,s,:))),'color',cmap_ppt(2,:),'LineWidth',1.5); hold on;
    xlim([kays(1) kays(end)])
    %ylim([-0.01 0.025])
    if (s==5)
        ylabel('Mean nu')
    end
    stamp('')
    title([loc ' Small'])
    
    %% Rec
    f3=figure(3);
    subplot(4,4,s)
    plot(kays,(squeeze(mrec(1,s,:))),'color',cmap_ppt(3,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mrec(2,s,:))),'color',cmap_ppt(1,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mrec(3,s,:))),'color',cmap_ppt(2,:),'LineWidth',1.5); hold on;
    xlim([kays(1) kays(end)])
    %ylim([-5 -1])
    if (s==5)
        ylabel('Mean recruitment')
    end
    stamp('')
    title([loc ' Small'])
    
    %% Con
    f4=figure(4);
    subplot(4,4,s)
    plot(kays,(squeeze(mcon(1,s,:))),'color',cmap_ppt(3,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mcon(2,s,:))),'color',cmap_ppt(1,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mcon(3,s,:))),'color',cmap_ppt(2,:),'LineWidth',1.5); hold on;
    xlim([kays(1) kays(end)])
    %ylim([0 0.025])
    if (s==5)
        ylabel('Mean consumption')
    end
    stamp('')
    title([loc ' Small'])
    
    %% Feeding level
    f5=figure(5);
    subplot(4,4,s)
    plot(kays,(squeeze(mlev(1,s,:))),'color',cmap_ppt(3,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mlev(2,s,:))),'color',cmap_ppt(1,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mlev(3,s,:))),'color',cmap_ppt(2,:),'LineWidth',1.5); hold on;
    xlim([kays(1) kays(end)])
    ylim([0 1])
    if (s==5)
        ylabel('Mean feeding level')
    end
    stamp('')
    title([loc ' Small'])
    
    %% Zoop frac con
    f6=figure(6);
    subplot(4,4,s)
    plot(kays,(squeeze(mzfrac(1,s,:))),'color',cmap_ppt(3,:),'LineWidth',1.5); hold on;
    plot(kays,(squeeze(mzfrac(2,s,:))),'color',cmap_ppt(1,:),'LineWidth',1.5); hold on;
    %     plot(kays,(squeeze(mzfrac(3,s,:))),'color',cmap_ppt(2,:),'LineWidth',1.5); hold on;
    xlim([kays(1) kays(end)])
    ylim([0 1.1])
    if (s==5)
        ylabel('Mean fraction consumed')
    end
    stamp('')
    title([loc ' Prey'])
    
end %spots

%%
print(f1,'-dpng',[fpath 'Climatol_ksat_125_2000_tot_mean_biomass_type_all_locs.png'])
print(f2,'-dpng',[fpath 'Climatol_ksat_125_2000_nu_small_all_locs.png'])
print(f3,'-dpng',[fpath 'Climatol_ksat_125_2000_rec_med_all_locs.png'])
print(f4,'-dpng',[fpath 'Climatol_ksat_125_2000_con_small_all_locs.png'])
print(f5,'-dpng',[fpath 'Climatol_ksat_125_2000_lev_small_all_locs.png'])
print(f6,'-dpng',[fpath 'Climatol_ksat_125_2000_frac_zoo_all_locs.png'])

%% Residence time
Fin  = NaN*ones(2,length(spots),np);
Pin  = NaN*ones(3,length(spots),np);
Din  = NaN*ones(3,length(spots),np);
Fout = NaN*ones(2,length(spots),np);
Pout = NaN*ones(3,length(spots),np);
Dout = NaN*ones(3,length(spots),np);
FB   = NaN*ones(2,length(spots),np);
PB   = NaN*ones(3,length(spots),np);
DB   = NaN*ones(3,length(spots),np);

for n = 1:length(kays)
    Ka = kays(n);
    tk = num2str(int64(100*Ka));
    close all
    sname = ['Climatol_ksat',tk,'_locs_'];
    cfile = char(dp);
    
    load([dpath sname harv '_lastyr_sum_mean_biom.mat']);
    
    %% Adjust size of things to be even 3 x locs
    % per biomass (rates)
    rFrec = Frec./Fmean;
    rPrec = Prec./Pmean;
    rDrec = Drec./Dmean;
    
    Frep = [zeros(1,length(spots)); Frep(1,:)];
    Prep = [zeros(2,length(spots)); Prep(1,:)];
    Drep = [zeros(2,length(spots)); Drep(1,:)];
    
    %% add gains and losses
    Fin(:,:,n) = max(0,Fmgr) + rFrec;
    Pin(:,:,n) = max(0,Pmgr) + rPrec;
    Din(:,:,n) = max(0,Dmgr) + rDrec;
    
    Fout(:,:,n) = Fgam + Frep + Fnat + Fpred + Ffrate;
    Pout(:,:,n) = Pgam + Prep + Pnat + Ppred + Pfrate;
    Dout(:,:,n) = Dgam + Drep + Dnat + Dpred + Dfrate;
    
    FB(:,:,n) = Fmean;
    PB(:,:,n) = Pmean;
    DB(:,:,n) = Dmean;
    
end

%%
Fres1 = 1 ./ Fin;
Pres1 = 1 ./ Pin;
Dres1 = 1 ./ Din;

Fres2 = 1 ./ Fout;
Pres2 = 1 ./ Pout;
Dres2 = 1 ./ Dout;

% Save values for all locs and all ksats that combo
save([dpath 'Climatol_locs_ksat_125_2000_res_time.mat'],'Fin','Fout','FB',...
    'Pin','Pout','PB','Din','Dout','DB');

%% Pcolor plots of residence times
jays = 0:3;
ays = [1 kays];
[agrid,jgrid]=meshgrid(ays,jays);

[nj,ni,nk] = size(PB);

res1F2 = NaN*ones(nj+1,16,nk+1);
res1P2 = NaN*ones(nj+1,16,nk+1);
res1D2 = NaN*ones(nj+1,16,nk+1);
res2F2 = NaN*ones(nj+1,16,nk+1);
res2P2 = NaN*ones(nj+1,16,nk+1);
res2D2 = NaN*ones(nj+1,16,nk+1);

res1F2(1:2,:,1:nk) =Fres1;
res1P2(1:nj,:,1:nk)=Pres1;
res1D2(1:nj,:,1:nk)=Dres1;
res2F2(1:2,:,1:nk) =Fres2;
res2P2(1:nj,:,1:nk)=Pres2;
res2D2(1:nj,:,1:nk)=Dres2;

% colors
cmBP=cmocean('speed');

%% Only use 3 domain examples: EBS(10), PUp(16), HOT(13)
sid = [10;16;13];

for s=1:length(sid)
    domain = sid(s);
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Res 1 F
    f2=figure(2);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(res1F2(:,domain,:)))
    colormap(cmBP)
    caxis([0 60])
    set(gca,'YTick',1:3,'YTickLabel',{'S','M','L'})
    if (s==2)
        title({loc; 'F res in (d)'})
    else
        title({loc; ''})
    end
    
    % Res 1 P
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(res1P2(:,domain,:)))
    colormap(cmBP)
    caxis([0 60])
    set(gca,'YTick',1:3,'YTickLabel',{'S','M','L'})
    if (s==2)
        title({loc; 'P res in (d)'})
    else
        title({loc; ''})
    end
    
    % Res 1 D
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(res1D2(:,domain,:)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([0 60])
    set(gca,'YTick',1:3,'YTickLabel',{'S','M','L'})
    if (s==2)
        title({loc; 'D res in (d)'})
    else
        title({loc; ''})
    end
    if (s==2)
        xlabel('ksat multiplier')
    end
    stamp('')
    
    %% Res 2 F
    f5=figure(5);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(res2F2(:,domain,:)))
    cmocean('speed')
    caxis([0 60])
    set(gca,'YTick',1:3,'YTickLabel',{'S','M','L'})
    if (s==2)
        title({loc; 'F res out (d)'})
    else
        title({loc; ''})
    end
    
    % Res 2 P
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(res2P2(:,domain,:)))
    cmocean('speed')
    caxis([0 60])
    set(gca,'YTick',1:3,'YTickLabel',{'S','M','L'})
    if (s==2)
        title({loc; 'P res out (d)'})
    else
        title({loc; ''})
    end
   
    % Res 2 D
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(res2D2(:,domain,:)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('speed')
    caxis([0 60])
    set(gca,'YTick',1:3,'YTickLabel',{'S','M','L'})
    if (s==2)
        title({loc; 'D res out (d)'})
    else
        title({loc; ''})
    end
    if (s==2)
        xlabel('ksat multiplier')
    end
    stamp('')
    
end %spots
%%
print(f2,'-dpng',[fpath 'Climatol_ksat_125_2000_Res1_3locs.png'])
print(f5,'-dpng',[fpath 'Climatol_ksat_125_2000_Res2_3locs.png'])

%% Line plots of in and out
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    f1=figure(1);
    subplot(4,4,s)
    plot(kays,squeeze(Fres1(1,s,:)),'r','LineWidth',2); hold on;
    plot(kays,squeeze(Fres2(1,s,:)),'--r','LineWidth',2); hold on;
    if (s==2)
        title({loc; 'SF res (d)'})
        legend('in','out')
    else
        title({loc; ''})
    end
    
    f3=figure(3);
    subplot(4,4,s)
    plot(kays,squeeze(Pres1(1,s,:)),'b','LineWidth',2); hold on;
    plot(kays,squeeze(Pres2(1,s,:)),'--b','LineWidth',2); hold on;
    if (s==2)
        title({loc; 'SP res (d)'})
        legend('in','out')
    else
        title({loc; ''})
    end
    
    f4=figure(4);
    subplot(4,4,s)
    plot(kays,squeeze(Dres1(1,s,:)),'k','LineWidth',2); hold on;
    plot(kays,squeeze(Dres2(1,s,:)),'--k','LineWidth',2); hold on;
    if (s==2)
        title({loc; 'SD res (d)'})
        legend('in','out')
    else
        title({loc; ''})
    end
end
% print(f1,'-dpng',[fpath 'Climatol_ksat_125_2000_SFres_all_locs.png'])
% print(f3,'-dpng',[fpath 'Climatol_ksat_125_2000_SPres_all_locs.png'])
% print(f4,'-dpng',[fpath 'Climatol_ksat_125_2000_SDres_all_locs.png'])

%% Line plots of in and out
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    f10=figure(10);
    subplot(4,4,s)
    plot(kays,squeeze(Fres1(1,s,:)),'r','LineWidth',2); hold on;
    plot(kays,squeeze(Fres2(1,s,:)),'--r','LineWidth',2); hold on;
    if (s==2)
        title({loc; 'SF res (d)'})
        legend('in','out')
    else
        title({loc; ''})
    end
    xlim([1 5])
    
    f30=figure(30);
    subplot(4,4,s)
    plot(kays,squeeze(Pres1(1,s,:)),'b','LineWidth',2); hold on;
    plot(kays,squeeze(Pres2(1,s,:)),'--b','LineWidth',2); hold on;
    if (s==2)
        title({loc; 'SP res (d)'})
        legend('in','out')
    else
        title({loc; ''})
    end
    xlim([1 5])
    
    f40=figure(40);
    subplot(4,4,s)
    plot(kays,squeeze(Dres1(1,s,:)),'k','LineWidth',2); hold on;
    plot(kays,squeeze(Dres2(1,s,:)),'--k','LineWidth',2); hold on;
    if (s==2)
        title({loc; 'SD res (d)'})
        legend('in','out')
    else
        title({loc; ''})
    end
    xlim([1 5])
end
print(f10,'-dpng',[fpath 'Climatol_ksat_125_500_SFres_all_locs.png'])
print(f30,'-dpng',[fpath 'Climatol_ksat_125_500_SPres_all_locs.png'])
print(f40,'-dpng',[fpath 'Climatol_ksat_125_500_SDres_all_locs.png'])


