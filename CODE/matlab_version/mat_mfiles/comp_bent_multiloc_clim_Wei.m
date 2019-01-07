% Calculate different skill metrics for global bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
datap = '/Volumes/GFDL/NC/Matlab_new_size/';
cp = '/Volumes/GFDL/CSV/Matlab_new_size/';

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

%% put on same grid as POEM output
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

test=seafl(:,2);
id=find(test>80);
test(id)=test(id)-360;

geolat_c = double(lat);
geolon_c = double(lon) - 280;

Zi = griddata(seafl(:,1),test,invert,geolat_c,geolon_c);
Zf = griddata(seafl(:,1),test,fish,geolat_c,geolon_c);

%% Maps
% Inv
figure(10)
surf(geolon_t,geolat_t,log10(Zi)); view(2); hold on;
shading flat
title('log10 mean benthic invert biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
%print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_inverts.png')

% Fish
figure(20)
surf(geolon_t,geolat_t,log10(Zf)); view(2); hold on;
shading flat
title('log10 mean benthic fish biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
%print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_fish.png')

%% Skill metrics

% Stowe et al 2009
% Univariate (model vs observ)
% 1. correlation coefficient [-1 1]; want close to 1
% 2. root mean square error; want close to 0; considering the magnitude rather than the direction
% 3. average error (bias); want close to 0; neg & pos discrepancies can cancel each other
% 4. average absolute error; want close to 0; considering the magnitude rather than the direction
% 5. modeling efficiency; want close to 1; <0 means avg obs better than model; how well a model predicts relative to avg obs

metrics={'r','RMSE','AE','AAE','MEF','Taylor nstd','tRMSD'};
metrics=metrics';

%Sometimes it is appropriate to log-transform the observations and predictions before
%calculating goodness-of-fit statistics so that differences between predicted and observed
%values will not be highly skewed and dominated by a small proportion of high values.

%
bees = 0.05:0.025:0.1;
cmaxs = 20:5:30;
iskill = NaN*ones(length(bees),length(cmaxs),7);
fskill = iskill;

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
        load([dpath 'Means_Climatol_All_fish03_' cfile '.mat'],...
            'md_mean','ld_mean','b_mean');
        
        %% Model bent & dem
        Zmd=NaN*ones(ni,nj);
        Zld=NaN*ones(ni,nj);
        Zb=NaN*ones(ni,nj);
        
        Zmd(ID)=md_mean;
        Zld(ID)=ld_mean;
        Zb(ID)=b_mean;
        
        Zd = Zmd+Zld;
        
        close all
        % Bent
        figure(1)
        surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
        shading flat
        title('log10 mean benthic biomass (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-2.5 0.5])
        stamp(cfile)
        %print('-dpng',[ppath 'Global_Bent.png'])
        
        % D
        figure(2)
        surf(geolon_t,geolat_t,log10(Zd)); view(2); hold on;
        shading flat
        title('log10 mean model M&L D biomass (g m^-^2)')
        colormap('jet')
        colorbar('h')
        caxis([-2.5 0.5])
        stamp(cfile)
        print('-dpng',[fpath 'Clim_All_fish03_global_Demersal.png'])
        
        
        %% Inverts
        obs = real(log(Zi(:)));
        model = real(log(Zb(:)));
        
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
        iskill(e,c,1) = num/den;
        
        % root mean square error
        num=nansum((p-o).^2);
        iskill(e,c,2) = sqrt(num/n);
        
        % average error
        iskill(e,c,3) = nansum(p-o) / n;
        
        % average absolute error
        iskill(e,c,4) = nansum(abs(p-o)) / n;
        
        % modeling efficiency
        num1=nansum((o-omean).^2);
        num2=nansum((p-o).^2);
        iskill(e,c,5) = (num1-num2)/num1;
        
        % Taylor normalized std
        iskill(e,c,6) = psig/osig;
        
        % total root mean square difference
        q1=nansum((o-p).^2);
        iskill(e,c,7) = sqrt(q1/n);
        
        
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
        fskill(e,c,1) = num/den;
        
        % root mean square error
        num=nansum((p-o).^2);
        fskill(e,c,2) = sqrt(num/n);
        
        % average error
        fskill(e,c,3) = nansum(p-o) / n;
        
        % average absolute error
        fskill(e,c,4) = nansum(abs(p-o)) / n;
        
        % modeling efficiency
        num1=nansum((o-omean).^2);
        num2=nansum((p-o).^2);
        fskill(e,c,5) = (num1-num2)/num1;
        
        % Taylor normalized std
        fskill(e,c,6) = psig/osig;
        
        % total root mean square difference
        q1=nansum((o-p).^2);
        fskill(e,c,7) = sqrt(q1/n);
        
    end
end
cfile2 = ['Dc_enc70-b200_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_noCC_RE00100'];
save([cp 'Bio_rates/' cfile2 '_global_skill_Wei.mat'],'invert',...
    'fish','metrics','iskill','fskill');

%% Plot results
nc = length(cmaxs);
ne = length(bees);
encs2 = [bees 0.125];
cmaxs2 = [cmaxs 35];
[cgrid,egrid]=meshgrid(cmaxs2,encs2);
i2  = NaN*ones(ne+1,nc+1,length(metrics));
i2(1:ne,1:nc,:) = iskill;
f2  = NaN*ones(ne+1,nc+1,length(metrics));
f2(1:ne,1:nc,:) = fskill;

cmap_ther = colormap(cmocean('thermal')); close all;
cmap_revt = flipud(cmap_ther);

for s=1:7
    f3=figure(3);
    subplot(3,3,s)
    pcolor(egrid,cgrid,squeeze(i2(:,:,s)))
    if (s==1 || s==5)
        cmocean('thermal')
    else
        colormap(cmap_revt)
    end
    colorbar
    %caxis([0 1])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',cmaxs,'YTickLabel',cmaxs)
    xlabel('BE')
    ylabel('Cmax coeff')
    if (s==2)
        title(['Inverts ' metrics{s}])
    else
        title(metrics{s})
    end
    
    f4=figure(4);
    subplot(3,3,s)
    pcolor(egrid,cgrid,squeeze(f2(:,:,s)))
    if (s==1 || s==5)
        cmocean('thermal')
    else
        colormap(cmap_revt)
    end
    colorbar
    %caxis([0 1])
    set(gca,'XTick',bees,'XTickLabel',bees,...
        'YTick',cmaxs,'YTickLabel',cmaxs)
    xlabel('BE')
    ylabel('Cmax coeff')
    if (s==2)
        title(['Fish ' metrics{s}])
    else
        title(metrics{s})
    end
end
print(f3,'-dpng',[figp cfile2 '_Clim_All_fish03_CmaxBE_Wei_invert_skill.png'])
print(f4,'-dpng',[figp cfile2 '_Clim_All_fish03_CmaxBE_Wei_fish_skill.png'])


