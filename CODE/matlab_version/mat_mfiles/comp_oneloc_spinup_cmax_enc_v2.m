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

%%
encs = linspace(10,100,10);
cmaxs = linspace(10,100,10);
nc = length(cmaxs);
ne = length(encs);
allF  = NaN*ones(length(spots),nc+1,ne+1);
allP  = NaN*ones(length(spots),nc+1,ne+1);
allD  = NaN*ones(length(spots),nc+1,ne+1);
mcon  = NaN*ones(3,nc+1,ne+1);
mlev  = NaN*ones(3,nc+1,ne+1);

cfile1 = ['Dc_enc-b200_m-b175-k08_fcrit20_c-b250',...
            '_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100_EncCmaxTests'];
load([datap 'Bio_rates/' cfile1 '.mat'],'allF','allP');
FPrat = allF./(allF+allP);
eepF(:,:,1) = allF(7,:,:);
pupF(:,:,1) = allF(11,:,:);
eepP(:,:,1) = allP(7,:,:);
pupP(:,:,1) = allP(11,:,:);
eepRat(:,:,1) = FPrat(7,:,:);
pupRat(:,:,1) = FPrat(11,:,:);
clear allF allP FPrat

cfile2 = ['Dc_enc-b200_m-b175-k09_fcrit20_c-b250',...
            '_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100_EncCmaxTests'];
load([datap 'Bio_rates/' cfile2 '.mat'],'allF','allP');
FPrat = allF./(allF+allP);
eepF(:,:,2) = allF(7,:,:);
pupF(:,:,2) = allF(11,:,:);
eepP(:,:,2) = allP(7,:,:);
pupP(:,:,2) = allP(11,:,:);
eepRat(:,:,2) = FPrat(7,:,:);
pupRat(:,:,2) = FPrat(11,:,:);
clear allF allP FPrat

cfile3 = ['Dc_enc-b200_m-b200-k09_fcrit20_c-b250',...
            '_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100_EncCmaxTests'];
load([datap 'Bio_rates/' cfile3 '.mat'],'allF','allP');
FPrat = allF./(allF+allP);
eepF(:,:,3) = allF(7,:,:);
pupF(:,:,3) = allF(11,:,:);
eepP(:,:,3) = allP(7,:,:);
pupP(:,:,3) = allP(11,:,:);
eepRat(:,:,3) = FPrat(7,:,:);
pupRat(:,:,3) = FPrat(11,:,:);
clear allF allP FPrat

fign = 'Dc_enc-b200_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100_EncCmaxTests';

%%
[egrid,cgrid]=meshgrid([cmaxs 110],[encs 110]);
kb={'b=0.175;k=0.08','b=0.175;k=0.09','b=0.2;k=0.09'};

for s=1:3
    %% Sum mean biom over stages
    f1=figure(1);
    subplot(3,3,3*s-2)
    pcolor(cgrid,egrid,squeeze(log10(eepF(:,:,s))))
    colorbar
    colormap('jet')
    caxis([-2 2])
    if (s==1)
        title({'log10 Mean F Biom'; kb{s}})
    else
        title(kb{s})
    end
    xlabel('Cmax coeff')
    ylabel('Enc coeff')
    
    subplot(3,3,3*s-1)
    pcolor(cgrid,egrid,squeeze(log10(eepP(:,:,s))))
    colorbar
    colormap('jet')
    caxis([-2 2])
    if (s==1)
        title({'log10 Mean P Biom'; 'EEP'})
    else
        title(kb{s})
    end
    xlabel('Cmax coeff')
    ylabel('Enc coeff')
    
    subplot(3,3,3*s)
    pcolor(cgrid,egrid,squeeze(eepRat(:,:,s)))
    colorbar
    colormap(cmap_color_rb)
    caxis([0 1])
    if (s==1)
        title({'Mean F:P'; kb{s}})
    else
        title(kb{s})
    end
    xlabel('Cmax coeff')
    ylabel('Enc coeff')
    
    
    f2=figure(2);
    subplot(3,3,3*s-2)
    pcolor(cgrid,egrid,squeeze(log10(pupF(:,:,s))))
    colorbar
    colormap('jet')
    caxis([-2 2])
    if (s==1)
        title({'log10 Mean F Biom'; kb{s}})
    else
        title(kb{s})
    end
    xlabel('Cmax coeff')
    ylabel('Enc coeff')
    
    subplot(3,3,3*s-1)
    pcolor(cgrid,egrid,squeeze(log10(pupP(:,:,s))))
    colorbar
    colormap('jet')
    caxis([-2 2])
    if (s==1)
        title({'log10 Mean P Biom'; 'PUp'})
    else
        title(kb{s})
    end
    xlabel('Cmax coeff')
    ylabel('Enc coeff')
    
    subplot(3,3,3*s)
    pcolor(cgrid,egrid,squeeze(pupRat(:,:,s)))
    colorbar
    colormap(cmap_color_rb)
    caxis([0 1])
    if (s==1)
        title({'Mean F:P'; kb{s}})
    else
        title(kb{s})
    end
    xlabel('Cmax coeff')
    ylabel('Enc coeff')
    %stamp(cfile2)
    
end %spots
print(f1,'-dpng',[figp sname fign '_FP_biomass_EEP.png'])
print(f2,'-dpng',[figp sname fign '_FP_biomass_PUp.png'])



