% Visualize output of POEM biol rate eq tests
% Spinup at one location
% 150 years, monthly means saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

efn = 70;
tefn = num2str(efn);
cfn = 20;
tcfn = num2str(cfn);
fcrit = 20;
nmort = '1';
kad = 50;
D = 'D075';
J = 'J100';
Ad = 'A050';
Sm = 'Sm025';
BE = '05';
rfrac = '00100';

kays = 0.0405:0.01:0.0916;
bees = 0.1:0.025:0.3;

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
harv = 'All_fish03';

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
load('MyColormaps.mat')
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
for C = 1:length(kays)
    kt = kays(C);
    for E = 1:length(bees)
        bpow = bees(E);
        tkfn = num2str(100+int64(100*kt));
        tbfn = num2str(1000+int64(1000*bpow));
        
        close all
%         dp = ['Dc_enc',tefn,'_cmax-metab',tcfn,'_b',tbfn(2:end),'_k',tkfn(2:end),...
%             '_fcrit',num2str(fcrit),'_',D,'_',J,'_',Ad,'_',Sm,...
%             '_nmort',nmort,'_BE',BE,'_CC',CC,'_lgRE',rfrac,'_mdRE',rfrac2];
        dp = ['Dc_enc70-b200_cm20_m-b',tbfn(2:end),'-k',tkfn(2:end),...
            '_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
        dpath = [datap char(dp) '/'];
        fpath = [figp char(dp) '/'];
        %         if (~isdir([figp char(dp)]))
        %             mkdir([figp char(dp)])
        %         end
        cfile = char(dp);
        
        all_mean=NaN*ones(3,4,length(spots));
        z = NaN*ones(length(spots),3);
        
        load([dpath sname 'locs_' harv '.mat'])
        SP = S_Sml_p;
        SF = S_Sml_f;
        SD = S_Sml_d;
        MP = S_Med_p;
        MF = S_Med_f;
        MD = S_Med_d;
        LP = S_Lrg_p;
        LD = S_Lrg_d;
        CO = S_Cobalt;
        
        t=1:size(SP,1);
        lyr=t((end-12+1):end);
        
        %% Final mean biomass in each size
        
        SP_mean=squeeze(mean(SP(lyr,1,:)))';
        SF_mean=squeeze(mean(SF(lyr,1,:)))';
        SD_mean=squeeze(mean(SD(lyr,1,:)))';
        MP_mean=squeeze(mean(MP(lyr,1,:)))';
        MF_mean=squeeze(mean(MF(lyr,1,:)))';
        MD_mean=squeeze(mean(MD(lyr,1,:)))';
        LP_mean=squeeze(mean(LP(lyr,1,:)))';
        LD_mean=squeeze(mean(LD(lyr,1,:)))';
        B_mean =squeeze(mean(CO(lyr,1,:)))';
        
        Pmean=[SP_mean;MP_mean;LP_mean];
        Fmean=[SF_mean;MF_mean];
        Dmean=[SD_mean;MD_mean;LD_mean];
        Bmean = B_mean;
        
        all_mean(1:2,1,:) = Fmean;
        all_mean(:,2,:) = Pmean;
        all_mean(:,3,:) = Dmean;
        all_mean(1,4,:) = Bmean;
        
        % Size spectrum (sum stages)
        spec = nansum(all_mean,2);
        
        %% Growth rate (nu - energy for biomass production)
        SP_mgr=squeeze(nanmean(SP(lyr,15,:)))';
        SF_mgr=squeeze(nanmean(SF(lyr,15,:)))';
        SD_mgr=squeeze(nanmean(SD(lyr,15,:)))';
        MP_mgr=squeeze(nanmean(MP(lyr,15,:)))';
        MF_mgr=squeeze(nanmean(MF(lyr,15,:)))';
        MD_mgr=squeeze(nanmean(MD(lyr,15,:)))';
        LP_mgr=squeeze(nanmean(LP(lyr,15,:)))';
        LD_mgr=squeeze(nanmean(LD(lyr,15,:)))';
        
        Pmgr=[SP_mgr;MP_mgr;LP_mgr];
        Fmgr=[SF_mgr;MF_mgr];
        Dmgr=[SD_mgr;MD_mgr;LD_mgr];
        
        %% Consump per biomass (I)
        SP_con=squeeze(nanmean(SP(lyr,14,:)))';
        SF_con=squeeze(nanmean(SF(lyr,14,:)))';
        SD_con=squeeze(nanmean(SD(lyr,14,:)))';
        MP_con=squeeze(nanmean(MP(lyr,14,:)))';
        MF_con=squeeze(nanmean(MF(lyr,14,:)))';
        MD_con=squeeze(nanmean(MD(lyr,14,:)))';
        LP_con=squeeze(nanmean(LP(lyr,14,:)))';
        LD_con=squeeze(nanmean(LD(lyr,14,:)))';
        
        Pcon=[SP_con;MP_con;LP_con];
        Fcon=[SF_con;MF_con];
        Dcon=[SD_con;MD_con;LD_con];
        
        %% Feeding level = con/cmax
        SP_lev=squeeze(nanmean(SP(lyr,20,:)))';
        SF_lev=squeeze(nanmean(SF(lyr,20,:)))';
        SD_lev=squeeze(nanmean(SD(lyr,20,:)))';
        MP_lev=squeeze(nanmean(MP(lyr,20,:)))';
        MF_lev=squeeze(nanmean(MF(lyr,20,:)))';
        MD_lev=squeeze(nanmean(MD(lyr,20,:)))';
        LP_lev=squeeze(nanmean(LP(lyr,20,:)))';
        LD_lev=squeeze(nanmean(LD(lyr,20,:)))';
        
        Plev=[SP_lev;MP_lev;LP_lev];
        Flev=[SF_lev;MF_lev];
        Dlev=[SD_lev;MD_lev;LD_lev];
        
        %% Fraction zoop losses consumed
        z(:,1) = squeeze(nanmean(CO(lyr,3,:)));
        z(:,2) = squeeze(nanmean(CO(lyr,4,:)));
        z(:,3) = squeeze(nanmean(CO(lyr,5,:)));
        
        %% Production (= nu * biom)
        SP_prod=squeeze(nanmean(SP(lyr,21,:)))';
        SF_prod=squeeze(nanmean(SF(lyr,21,:)))';
        SD_prod=squeeze(nanmean(SD(lyr,21,:)))';
        MP_prod=squeeze(nanmean(MP(lyr,21,:)))';
        MF_prod=squeeze(nanmean(MF(lyr,21,:)))';
        MD_prod=squeeze(nanmean(MD(lyr,21,:)))';
        LP_prod=squeeze(nanmean(LP(lyr,21,:)))';
        LD_prod=squeeze(nanmean(LD(lyr,21,:)))';
        
        Pprod=[SP_prod;MP_prod;LP_prod];
        Fprod=[SF_prod;MF_prod];
        Dprod=[SD_prod;MD_prod;LD_prod];
        
        %% Reproduction
        Frep(1,:)=squeeze(nanmean(MF(lyr,18,:)))';
        Drep(1,:)=squeeze(nanmean(LD(lyr,18,:)))';
        Prep(1,:)=squeeze(nanmean(LP(lyr,18,:)))';
        Frep(2,:)=squeeze(nanmean(MF(lyr,1,:).*MF(lyr,18,:)))';
        Drep(2,:)=squeeze(nanmean(LD(lyr,1,:).*LD(lyr,18,:)))';
        Prep(2,:)=squeeze(nanmean(LP(lyr,1,:).*LP(lyr,18,:)))';
        
        %% Metabolism
        SP_met=squeeze(nanmean(SP(lyr,24,:)))';
        SF_met=squeeze(nanmean(SF(lyr,24,:)))';
        SD_met=squeeze(nanmean(SD(lyr,24,:)))';
        MP_met=squeeze(nanmean(MP(lyr,24,:)))';
        MF_met=squeeze(nanmean(MF(lyr,24,:)))';
        MD_met=squeeze(nanmean(MD(lyr,24,:)))';
        LP_met=squeeze(nanmean(LP(lyr,24,:)))';
        LD_met=squeeze(nanmean(LD(lyr,24,:)))';
        
        Pmet=[SP_met;MP_met;LP_met];
        Fmet=[SF_met;MF_met];
        Dmet=[SD_met;MD_met;LD_met];
        
        %% Predation
        SP_pred=squeeze(nanmean(SP(lyr,22,:)))';
        SF_pred=squeeze(nanmean(SF(lyr,22,:)))';
        SD_pred=squeeze(nanmean(SD(lyr,22,:)))';
        MP_pred=squeeze(nanmean(MP(lyr,22,:)))';
        MF_pred=squeeze(nanmean(MF(lyr,22,:)))';
        MD_pred=squeeze(nanmean(MD(lyr,22,:)))';
        LP_pred=squeeze(nanmean(LP(lyr,22,:)))';
        LD_pred=squeeze(nanmean(LD(lyr,22,:)))';
        
        Ppred=[SP_pred;MP_pred;LP_pred];
        Fpred=[SF_pred;MF_pred];
        Dpred=[SD_pred;MD_pred;LD_pred];
        
        %% Natural mortality
        Pnat(1,:)=squeeze(nanmean(SP(lyr,23,:)))';
        Fnat(1,:)=squeeze(nanmean(SF(lyr,23,:)))';
        Dnat(1,:)=squeeze(nanmean(SD(lyr,23,:)))';
        Pnat(2,:)=squeeze(nanmean(MP(lyr,23,:)))';
        Fnat(2,:)=squeeze(nanmean(MF(lyr,23,:)))';
        Dnat(2,:)=squeeze(nanmean(MD(lyr,23,:)))';
        Pnat(3,:)=squeeze(nanmean(LP(lyr,23,:)))';
        Dnat(3,:)=squeeze(nanmean(LD(lyr,23,:)))';
        
        %% Fishing
        MP_fish=squeeze(nanmean(MP(lyr,25,:)))';
        MF_fish=squeeze(nanmean(MF(lyr,25,:)))';
        MD_fish=squeeze(nanmean(MD(lyr,25,:)))';
        LP_fish=squeeze(nanmean(LP(lyr,25,:)))';
        LD_fish=squeeze(nanmean(LD(lyr,25,:)))';
        
        Pfish=[zeros(size(MP_fish));MP_fish;LP_fish];
        Ffish=[zeros(size(MF_fish));MF_fish];
        Dfish=[zeros(size(MD_fish));MD_fish;LD_fish];
        
        MP_totcatch=squeeze(nansum(MP(lyr,25,:)))';
        MF_totcatch=squeeze(nansum(MF(lyr,25,:)))';
        MD_totcatch=squeeze(nansum(MD(lyr,25,:)))';
        LP_totcatch=squeeze(nansum(LP(lyr,25,:)))';
        LD_totcatch=squeeze(nansum(LD(lyr,25,:)))';
        
        Ptotcatch=[zeros(size(MP_totcatch));MP_totcatch;LP_totcatch];
        Ftotcatch=[zeros(size(MF_totcatch));MF_totcatch];
        Dtotcatch=[zeros(size(MD_totcatch));MD_totcatch;LD_totcatch];
        
        %% Total mortality w/o fishing
        Fmort = Fpred + Fnat;
        Pmort = Ppred + Pnat;
        Dmort = Dpred + Dnat;
        
        %% Total mortality w/ fishing
        Fmortf = Fpred + Fnat + Ffish;
        Pmortf = Ppred + Pnat + Pfish;
        Dmortf = Dpred + Dnat + Dfish;
        
        %% Gross growth efficiency (= nu/consump)
        SP_gge=squeeze(nanmean(SP(lyr,15,:)./SP(lyr,14,:)))';
        SF_gge=squeeze(nanmean(SF(lyr,15,:)./SF(lyr,14,:)))';
        SD_gge=squeeze(nanmean(SD(lyr,15,:)./SD(lyr,14,:)))';
        MP_gge=squeeze(nanmean(MP(lyr,15,:)./MP(lyr,14,:)))';
        MF_gge=squeeze(nanmean(MF(lyr,15,:)./MF(lyr,14,:)))';
        MD_gge=squeeze(nanmean(MD(lyr,15,:)./MD(lyr,14,:)))';
        LP_gge=squeeze(nanmean(LP(lyr,15,:)./LP(lyr,14,:)))';
        LD_gge=squeeze(nanmean(LD(lyr,15,:)./LD(lyr,14,:)))';
        
        Pgge=[SP_gge;MP_gge;LP_gge];
        Fgge=[SF_gge;MF_gge];
        Dgge=[SD_gge;MD_gge;LD_gge];
        
        save([dpath sname '_' harv '_lastyr_sum_mean_biom.mat'],'Pmean','Fmean','Dmean','all_mean',...
            'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
            'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
            'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
            'Dtotcatch','Pgge','Fgge','Dgge','Plev','Flev','Dlev','Bmean');
        
    end %enc
end %cmax

%%
nc = length(kays);
ne = length(bees);
allF  = NaN*ones(length(spots),nc+1,ne+1);
allP  = NaN*ones(length(spots),nc+1,ne+1);
allD  = NaN*ones(length(spots),nc+1,ne+1);
mcon  = NaN*ones(3,nc+1,ne+1);
mlev  = NaN*ones(3,nc+1,ne+1);
for C = 1:length(kays)
    kt = kays(C);
    for E = 1:length(bees)
        bpow = bees(E);
        tkfn = num2str(100+int64(100*kt));
        tbfn = num2str(1000+int64(1000*bpow));
        
        close all
        %         dp = ['Dc_enc',tefn,'_cmax-metab',tcfn,'_b',tbfn(2:end),'_k',tkfn(2:end),...
%             '_fcrit',num2str(fcrit),'_',D,'_',J,'_',Ad,'_',Sm,...
%             '_nmort',nmort,'_BE',BE,'_CC',CC,'_lgRE',rfrac,'_mdRE',rfrac2];
        dp = ['Dc_enc70-b200_cm20_m-b',tbfn(2:end),'-k',tkfn(2:end),...
            '_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
        dpath = [datap char(dp) '/'];
        fpath = [figp char(dp) '/'];
        cfile = char(dp);
        cfile2 = ['Dc_enc',tefn,'_cmax-metab',tcfn,...
            '_fcrit',num2str(fcrit),'_',D,'_',J,'_',Ad,'_',Sm,...
            '_nmort',nmort,'_BE',BE,'_noCC_RE',rfrac,...
            '_CmaxBKtests'];
        
        load([dpath sname '_' harv '_lastyr_sum_mean_biom.mat']);
        
        %% Biomass of each type
        allF(:,C,E) = squeeze(nansum(all_mean(:,1,:)));
        allP(:,C,E) = squeeze(nansum(all_mean(:,2,:)));
        allD(:,C,E) = squeeze(nansum(all_mean(:,3,:)));
        
        %% Consump vs. weight
        Pcon = Pcon .* mass .* 365;
        Fcon = Fcon .* mass(1:2,:) .* 365;
        Dcon = Dcon .* mass .* 365;
        Fcon(3,:) = nan;
        Tcon = [Fcon,Pcon,Dcon];
        Mcon = nanmean(Tcon,2);
        mcon(:,C,E) = Mcon;
        
        %% Feeding level
        Flev(3,:) = nan;
        Tlev = [Flev,Plev,Dlev];
        Mlev = nanmean(Tlev,2);
        mlev(:,C,E) = Mlev;
        
    end %enc
end %kays

%% Save values for all locs and all kays that combo
save([datap 'Bio_rates/' cfile2 '.mat'],'allF','allP','allD','mcon','mlev');

%%
kays2 = [kays 0.105];
bees2 = [bees 0.325];
[bgrid,kgrid]=meshgrid(bees2,kays2);
FPrat = allF./(allF+allP);
DPrat = allD./(allD+allP);
for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f2=figure(2);
    subplot(4,3,s)
    pcolor(bgrid,kgrid,squeeze(log10(allF(s,:,:))))
    colorbar
    caxis([-2 2])
    set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
    if (s==2)
        title({'log10 Mean F Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('Weight exponent')
    ylabel('Temp exponent')
    stamp(cfile2)
    
    f3=figure(3);
    subplot(4,3,s)
    pcolor(bgrid,kgrid,squeeze(log10(allP(s,:,:))))
    colorbar
    caxis([-2 2])
    set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
    if (s==2)
        title({'log10 Mean P Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('Weight exponent')
    ylabel('Temp exponent')
    stamp(cfile2)
    
    
    f4=figure(4);
    subplot(4,3,s)
    pcolor(bgrid,kgrid,squeeze(log10(allD(s,:,:))))
    colorbar
    caxis([-2 2])
    set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
    if (s==2)
        title({'log10 Mean D Biom (g m^-^2) in final year'; loc})
    else
        title(loc)
    end
    xlabel('Weight exponent')
    ylabel('Temp exponent')
    stamp(cfile2)
    
    f5=figure(5);
    subplot(4,3,s)
    pcolor(bgrid,kgrid,squeeze(FPrat(s,:,:)))
    colorbar
    colormap(cmap_color_rb)
    caxis([0 1])
    set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
    if (s==2)
        title({'Frac F:P in final year'; loc})
    else
        title(loc)
    end
    xlabel('Weight exponent')
    ylabel('Temp exponent')
    stamp(cfile2)
    
    f6=figure(6);
    subplot(4,3,s)
    pcolor(bgrid,kgrid,squeeze(DPrat(s,:,:)))
    colorbar
    colormap(cmap_color_rb)
    caxis([0 1])
    set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
    if (s==2)
        title({'Frac D:P in final year'; loc})
    else
        title(loc)
    end
    xlabel('Weight exponent')
    ylabel('Temp exponent')
    stamp(cfile2)
    
end %spots
print(f2,'-dpng',[figp sname cfile2 '_totF_mean_biomass_type_all_locs.png'])
print(f3,'-dpng',[figp sname cfile2 '_totP_mean_biomass_type_all_locs.png'])
print(f4,'-dpng',[figp sname cfile2 '_totD_mean_biomass_type_all_locs.png'])
print(f5,'-dpng',[figp sname cfile2 '_FP_frac_all_locs.png'])
print(f6,'-dpng',[figp sname cfile2 '_DP_frac_all_locs.png'])

%% Feeding level
figure(10)
subplot(2,2,1)
pcolor(bgrid,kgrid,squeeze(mlev(1,:,:)))
colorbar
caxis([0 1])
set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
xlabel('Weight exponent')
ylabel('Temp exponent')
title('S Mean feeding level')

subplot(2,2,2)
pcolor(bgrid,kgrid,squeeze(mlev(2,:,:)))
colorbar
caxis([0 1])
set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
xlabel('Weight exponent')
ylabel('Temp exponent')
title('M Mean feeding level')

subplot(2,2,3)
pcolor(bgrid,kgrid,squeeze(mlev(3,:,:)))
colorbar
caxis([0 1])
set(gca,'XTick',bees(2:2:end),'XTickLabel',bees(1:2:end),...
        'YTick',kays(2:2:end),'YTickLabel',kays(1:2:end))
xlabel('Weight exponent')
ylabel('Temp exponent')
title('L Mean feeding level')
print('-dpng',[figp sname cfile2 '_mean_flev_size_all_locs.png'])
