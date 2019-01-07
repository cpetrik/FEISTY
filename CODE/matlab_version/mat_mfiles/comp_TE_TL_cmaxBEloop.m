% Compare TEs and TLs

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';
coast = [3;4;6;8;10];
coastal = spots(coast)';


% POEM file info
frate = 0.3;
tfish = num2str(100+int64(10*frate));
sname = 'Clim_';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

stages={'SF','MF','SP','MP','LP','SD','MD','LD'};

%%
bees = 0.05:0.025:0.1;
cmaxs = 20:5:30;
TL = NaN*ones(length(bees),length(cmaxs),length(spots));
TEc = TL;
TEp1 = TL;
TEp2 = TL;
TEe1 = TL;
TEe2 = TL;
for e=1:length(bees)
    for c=1:length(cmaxs)
        bent_eff = bees(e);
        h = cmaxs(c);
        tcfn = num2str(h);
        tbe = num2str(100+int64(100*bent_eff));
        
        cfile = ['Dc_enc70-b200_cm',tcfn,...
            '_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE',...
            tbe(2:end),'_noCC_RE00100'];
        dpath = [datap cfile '/'];
        fpath = [figp cfile '/'];
        load([dpath sname 'locs_' harv '_lastyr_TE_TL.mat'])
        
        %%
        TL(e,c,:) = TLprey(8,:);        %Dem
        TEc(e,c,:) = TEcon(2,:);        %Large
        TEp1(e,c,:) = TEprod_grp(3,:);  %Large
        TEp2(e,c,:) = TEprod(5,:);      %Dem
        TEe1(e,c,:) = TEeff1(2,:);      %Large
        TEe2(e,c,:) = TEeff2(2,:);      %Large
            
    end
end
cfile2 = ['Dc_enc70-b200_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_noCC_RE00100'];
save([datap 'Bio_rates/' cfile2 '_locs_CmaxBEtests_TE_TL.mat']);

%% Plots
nc = length(cmaxs);
ne = length(bees);
encs2 = [bees 0.125];
cmaxs2 = [cmaxs 35];
[cgrid,egrid]=meshgrid(cmaxs2,encs2);
r2  = NaN*ones(ne+1,nc+1,length(spots));
r2(1:ne,1:nc,:) = TL;
te2  = NaN*ones(ne+1,nc+1,length(spots));
te2(1:ne,1:nc,:) = TEp2;

cmap_ther = colormap(cmocean('thermal')); close all;
cmap_revt = flipud(cmap_ther);
for s=1:length(spots)
    f1=figure(1);
    subplot(4,4,s)
    pcolor(egrid,cgrid,squeeze(r2(:,:,s)))
    cmocean('thermal')
    colorbar
    caxis([3 3.5])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',cmaxs,'YTickLabel',cmaxs)
    xlabel('BE')
    ylabel('Cmax coeff')
    if (s~=2)
        title(spots{s})
    else
        title([spots{s} ' LD TL'])
    end
    
    f2=figure(2);
    subplot(4,4,s)
    pcolor(egrid,cgrid,squeeze(te2(:,:,s)))
    cmocean('thermal')
    colorbar
    caxis([0 0.005])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',cmaxs,'YTickLabel',cmaxs)
    xlabel('BE')
    ylabel('Cmax coeff')
    if (s~=2)
        title(spots{s})
    else
        title([spots{s} ' LD TE'])
    end
end
print(f1,'-dpng',[figp cfile2 '_CmaxBE_LD_TL.png'])
print(f2,'-dpng',[figp cfile2 '_CmaxBE_LD_TEprod.png'])
