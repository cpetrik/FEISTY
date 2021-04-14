% Calc residence time
% Residence time = 1 / input
% or             = 1 / output
% Total inputs: rec, nu
% Total outputs: gamma, rep, nmort, die (pred), yield (fishing)
% Use per biomass rates (g/g/m2/d)

clear all
close all

%%
cfile = 'Dc_enc-b200_m4-b175-k086_c-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/MIP/NC/Matlab_new_size/param_ensemble/',cfile,'/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/param_ensemble/';
ppath = [pp cfile '/'];

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught','fmort'};
cols=cols';
spots=spots';
spots{16} = 'PUP';

%%
load([fpath 'Locs_Climatol_All_fish03_means_aenc_acmax_search.mat']);

%% add gains and losses
% +nu is g/g/m2/d; rec is g/m2/d
SF_in = max(0,nSF) + rSF ./ SF;
SP_in = max(0,nSP) + rSP ./ SP;
SD_in = max(0,nSD) + rSD ./ SD;
MF_in = max(0,nMF) + rMF ./ MF;
MP_in = max(0,nMP) + rMP ./ MP;
MD_in = max(0,nMD) + rMD ./ MD;
LP_in = max(0,nLP) + rLP ./ LP;
LD_in = max(0,nLD) + rLD ./ LD;

% (-nu,) gamma, rep, nmort, pred, fmort are g/g/m2/d; mort, die, yieLD are g/m2/d
SF_out = abs(min(0,nSF)) + gSF + mSF + pSF;
SP_out = abs(min(0,nSP)) + gSP + mSP + pSP;
SD_out = abs(min(0,nSD)) + gSD + mSD + pSD;
MF_out = abs(min(0,nMF)) + gMF + mMF + pMF + fMF + eMF;
MP_out = abs(min(0,nMP)) + gMP + mMP + pMP + fMP;
MD_out = abs(min(0,nMD)) + gMD + mMD + pMD + fMD;
LP_out = abs(min(0,nLP)) + gLP + mLP + pLP + fLP + eLP;
LD_out = abs(min(0,nLD)) + gLD + mLD + pLD + fLD + eLD;

%% v1
SF_res1 = 1 ./ SF_in;
SP_res1 = 1 ./ SP_in;
SD_res1 = 1 ./ SD_in;
MF_res1 = 1 ./ MF_in;
MP_res1 = 1 ./ MP_in;
MD_res1 = 1 ./ MD_in;
LP_res1 = 1 ./ LP_in;
LD_res1 = 1 ./ LD_in;

% SF_res1(isinf(SF_res1(:))) = NaN;
% SP_res1(isinf(SP_res1(:))) = NaN;
% SD_res1(isinf(SD_res1(:))) = NaN;
% MF_res1(isinf(MF_res1(:))) = NaN;
% MP_res1(isinf(MP_res1(:))) = NaN;
% MD_res1(isinf(MD_res1(:))) = NaN;
% LP_res1(isinf(LP_res1(:))) = NaN;
% LD_res1(isinf(LD_res1(:))) = NaN;

%% v2
SF_res2 = 1 ./ SF_out;
SP_res2 = 1 ./ SP_out;
SD_res2 = 1 ./ SD_out;
MF_res2 = 1 ./ MF_out;
MP_res2 = 1 ./ MP_out;
MD_res2 = 1 ./ MD_out;
LP_res2 = 1 ./ LP_out;
LD_res2 = 1 ./ LD_out;

%% Save
save([fpath 'Residence_time_means_Climatol_' harv '_locs_aenc_acmax.mat'],...
  'SF_res1','SP_res1','SD_res1','MF_res1','MP_res1','MD_res1','LP_res1','LD_res1',...
  'SF_res2','SP_res2','SD_res2','MF_res2','MP_res2','MD_res2','LP_res2','LD_res2')

%% Histograms
edges = [0:30:360 547 730 912 1095];
edge2 = -4:2;

%biomass
figure(10)
subplot(3,3,1)
histogram(log10(SF),edge2)
title('SF')

subplot(3,3,2)
histogram(log10(SP),edge2)
title('SP')

subplot(3,3,3)
histogram(log10(SD),edge2)
title('SD')

subplot(3,3,4)
histogram(log10(MF),edge2)
title('MF')

subplot(3,3,5)
histogram(log10(MP),edge2)
title('MP')

subplot(3,3,6)
histogram(log10(MD),edge2)
title('MD')

subplot(3,3,8)
histogram(log10(LP),edge2)
title('LP')

subplot(3,3,9)
histogram(log10(LD),edge2)
title('LD')

%res 1
figure(11)
subplot(3,3,1)
histogram((SF_res1),edges)
title('SF')

subplot(3,3,2)
histogram((SP_res1),edges)
title('SP')

subplot(3,3,3)
histogram((SD_res1),edges)
title('SD')

subplot(3,3,4)
histogram((MF_res1),edges)
title('MF')

subplot(3,3,5)
histogram((MP_res1),edges)
title('MP')

subplot(3,3,6)
histogram((MD_res1),edges)
title('MD')

subplot(3,3,8)
histogram((LP_res1),edges)
title('LP')

subplot(3,3,9)
histogram((LD_res1),edges)
title('LD')

%res 2
figure(12)
subplot(3,3,1)
histogram((SF_res2),edges)
title('SF')

subplot(3,3,2)
histogram((SP_res2),edges)
title('SP')

subplot(3,3,3)
histogram((SD_res2),edges)
title('SD')

subplot(3,3,4)
histogram((MF_res2),edges)
title('MF')

subplot(3,3,5)
histogram((MP_res2),edges)
title('MP')

subplot(3,3,6)
histogram((MD_res2),edges)
title('MD')

subplot(3,3,8)
histogram((LP_res2),edges)
title('LP')

subplot(3,3,9)
histogram((LD_res2),edges)
title('LD')

%% Figures
nj=length(aep);
nk=length(acp);

jays = [aep 70];
ays = [acp 20];
[agrid,jgrid]=meshgrid(ays,jays);

mSFres2 = NaN*ones(nj+1,nk+1,16);
mSPres2 = NaN*ones(nj+1,nk+1,16);
mSDres2 = NaN*ones(nj+1,nk+1,16);
mMFres2 = NaN*ones(nj+1,nk+1,16);
mMPres2 = NaN*ones(nj+1,nk+1,16);
mMDres2 = NaN*ones(nj+1,nk+1,16);
mLPres2 = NaN*ones(nj+1,nk+1,16);
mLDres2 = NaN*ones(nj+1,nk+1,16);
mSFbio = NaN*ones(nj+1,nk+1,16);
mSPbio = NaN*ones(nj+1,nk+1,16);
mSDbio = NaN*ones(nj+1,nk+1,16);
mMFbio = NaN*ones(nj+1,nk+1,16);
mMPbio = NaN*ones(nj+1,nk+1,16);
mMDbio = NaN*ones(nj+1,nk+1,16);
mLPbio = NaN*ones(nj+1,nk+1,16);
mLDbio = NaN*ones(nj+1,nk+1,16);

mSFres2(1:nj,1:nk,:)=SF_res2;
mSPres2(1:nj,1:nk,:)=SP_res2;
mSDres2(1:nj,1:nk,:)=SD_res2;
mMFres2(1:nj,1:nk,:)=MF_res2;
mMPres2(1:nj,1:nk,:)=MP_res2;
mMDres2(1:nj,1:nk,:)=MD_res2;
mLPres2(1:nj,1:nk,:)=LP_res2;
mLDres2(1:nj,1:nk,:)=LD_res2;
mSFbio(1:nj,1:nk,:)=SF;
mSPbio(1:nj,1:nk,:)=SP;
mSDbio(1:nj,1:nk,:)=SD;
mMFbio(1:nj,1:nk,:)=MF;
mMPbio(1:nj,1:nk,:)=MP;
mMDbio(1:nj,1:nk,:)=MD;
mLPbio(1:nj,1:nk,:)=LP;
mLDbio(1:nj,1:nk,:)=LD;

Fbio = mSFbio + mMFbio ;
Pbio = mSPbio + mMPbio + mLPbio;
Dbio = mSDbio + mMDbio + mLDbio;

%%

%% Only use 3 domain examples: EBS(10), PUp(16), HOT(13)
sid = [10;16;13];
for s=1:length(sid)
    domain = sid(s);
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Biomass
    f1=figure(1);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(log10(Fbio(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('matter')
    caxis([-2 2])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'log10 Mean F Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(log10(Pbio(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('matter')
    caxis([-2 2])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('log10 Mean P Biom (g m^-^2)')
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(log10(Dbio(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('matter')
    caxis([-2 2])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('log10 Mean D Biom (g m^-^2)')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    %stamp(pname)
    
    %% Small res
    f2=figure(2);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(mSFres2(:,:,domain)))
    cmocean('speed')
    caxis([1 60])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'SF res'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(mSPres2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('speed')
    caxis([1 60])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('SP res')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(mSDres2(:,:,domain)))
    cmocean('speed')
    caxis([1 60])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('SD res')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    
    %% Medium res
    f3=figure(3);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(mMFres2(:,:,domain)))
    cmocean('speed')
    caxis([60 350])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'MF res'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(mMPres2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('speed')
    caxis([60 350])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('MP res')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(mMDres2(:,:,domain)))
    cmocean('speed')
    caxis([60 350])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('MD res')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    
    %% Large res
    f4=figure(4);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(mLPres2(:,:,domain)))
    cmocean('speed')
    caxis([90 600])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'LP res'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(mLDres2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('speed')
    caxis([90 600])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('LD res')
    end
    xlabel('a_C')
    if (s==1)
        ylabel('a_E')
    end
    
    
end %spots
%%
print(f1,'-dpng',[ppath 'biomass_locs_aenc_acmax_3locs.png'])
print(f2,'-dpng',[ppath 'Sres2_locs_aenc_acmax_3locs.png'])
print(f3,'-dpng',[ppath 'Mres2_locs_aenc_acmax_3locs.png'])
print(f4,'-dpng',[ppath 'Lres2_locs_aenc_acmax_3locs.png'])

%% check biomass in all locations
for s=1:length(spots)
    domain = s;
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Biomass
    f5=figure(5);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(Fbio(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('matter')
    caxis([-2 2])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'log10 Mean F Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    f6=figure(6);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(Pbio(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('matter')
    caxis([-2 2])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'log10 Mean P Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    f7=figure(7);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(log10(Dbio(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    cmocean('matter')
    caxis([-2 2])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'log10 Mean D Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    %% Small res
    f8=figure(8);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(mSFres2(:,:,domain)))
    cmocean('speed')
    caxis([1 60])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'SF res'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
    f9=figure(9);
    subplot(4,4,s)
    pcolor(agrid,jgrid,squeeze(mSPres2(:,:,domain)))
    cmocean('speed')
    caxis([1 60])
    set(gca,'XTick',(1:2:15)+5,'XTickLabel',5:2:20,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'SP res'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_E')
    end
    
end %spots
%%
print(f5,'-dpng',[ppath 'Fbiomass_locs_aenc_acmax_16locs.png'])
print(f6,'-dpng',[ppath 'Pbiomass_locs_aenc_acmax_16locs.png'])
print(f7,'-dpng',[ppath 'Dbiomass_locs_aenc_acmax_16locs.png'])
print(f8,'-dpng',[ppath 'SFres_locs_aenc_acmax_16locs.png'])
print(f9,'-dpng',[ppath 'SPres_locs_aenc_acmax_16locs.png'])



