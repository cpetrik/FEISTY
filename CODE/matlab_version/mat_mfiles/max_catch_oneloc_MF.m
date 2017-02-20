% Find RE of max LD catch for each loc and mort

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_big_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Fishing/';

fsim = {'.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
frate = {'01','02','03','04','05','06','07','08','09','10'};
Frate = 1:10;
RE = {'1000','0500','0100'};
reff = [1.0,0.5,0.1];
sreff = {'1.0','0.5','0.1'};
nmrt = [0,2:5];
Mort = {'None','Hartvig','mizer','J&C','P&W'};
fcrit = 40;
kad = 100;
pref = 'D100';
BE = '05';
CC = '050';

cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_BE' BE '_CC' CC '_MF_fishing_catch'];

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

%% colors
cm7=[0 0 0.75;... %b
    1 0 0;...     %r
    0 0 0;...     %black
    0 0.7 0;...   %g
    1 0 1;...     %m
    0/255 206/255 209/255;... %turq
    0.5 0 0];   %maroon     
    
set(groot,'defaultAxesColorOrder',cm7);

%%
maxC = NaN*ones(length(frate),length(RE),length(spots));
BF = NaN*ones(length(frate),length(RE),length(spots));
B0 = NaN*ones(length(RE),length(spots),length(nmrt));
for n = 1:length(nmrt)
    nmort = num2str(nmrt(n));
    for r = 1:length(RE)
        
        rfrac = RE{r};
        
        dp0 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac];
        dpath0 = [char(dp0) '/'];
        load([datap dpath0 sname 'lastyr_sum_mean_biom'],'all_mean');
        B0(r,:,n) = squeeze(all_mean(2,1,:));
        clear all_mean
        
        %%
        for i=1:length(frate)
            F = frate{i};
            dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_MF_fish' F];
            dpath = [char(dp) '/'];
            load([datap dpath sname 'lastyr_sum_mean_biom'],'Ftotcatch','all_mean');
            
            %% Fishing
            Ftot = sum(Ftotcatch);
            maxC(i,r,:) = Ftot;
            
            %% Unharvested Biomass
            BF(i,r,:) = squeeze(all_mean(2,1,:));
            
        end %frate
    end %RE
    
    if (n==1)
        m0_maxC = maxC;
        m0_BF = BF;
    elseif (n==2)
        m2_maxC = maxC;
        m2_BF = BF;
    elseif (n==3)
        m3_maxC = maxC;
        m3_BF = BF;
    elseif (n==4)
        m4_maxC = maxC;
        m4_BF = BF;
    elseif (n==5)
        m5_maxC = maxC;
        m5_BF = BF;
    end
end %nmort

%% Plot F catch
leg=1:5;
leg=repmat(leg,10,1);

for s = 1:length(spots)
    f1=figure(1);
    subplot(4,3,s)
    plot(1:10,m0_maxC(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{1}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total F catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print(f1,'-dpng',[figp cfile2 '_allRE_nmort0.png'])

for s = 1:length(spots)
    f2=figure(2);
    subplot(4,3,s)
    plot(1:10,m2_maxC(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{2}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total F catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print(f2,'-dpng',[figp cfile2 '_allRE_nmort2.png'])

for s = 1:length(spots)
    f3=figure(3);
    subplot(4,3,s)
    plot(1:10,m3_maxC(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{3}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total F catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print(f3,'-dpng',[figp cfile2 '_allRE_nmort3.png'])

for s = 1:length(spots)
    f4=figure(4);
    subplot(4,3,s)
    plot(1:10,m4_maxC(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{4}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total F catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print(f4,'-dpng',[figp cfile2 '_allRE_nmort4.png'])

for s = 1:length(spots)
    f5=figure(5);
    subplot(4,3,s)
    plot(1:10,m5_maxC(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{5}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total F catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print(f5,'-dpng',[figp cfile2 '_allRE_nmort5.png'])

%% Plot Bf:B0
cfile3 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_BE' BE '_CC' CC '_MF_BfB0'];
y = 0.4*ones(10,1);            

b0 = NaN*ones(length(frate),length(RE),length(spots));
for f=1:length(frate)
    b0(f,:,:) = B0(:,:,1);
end
bfb0 = m0_BF ./ b0;
for s = 1:length(spots)
    figure(6);
    subplot(4,3,s)
    plot(1:10,y,'--k'); hold on;
    plot(1:10,bfb0(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{1}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF B:B0 in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print('-dpng',[figp cfile3 '_allRE_nmort0.png'])

b0 = NaN*ones(length(frate),length(RE),length(spots));
for f=1:length(frate)
    b0(f,:,:) = B0(:,:,2);
end
bfb0 = m2_BF ./ b0;
for s = 1:length(spots)
    figure(7);
    subplot(4,3,s)
    plot(1:10,y,'--k'); hold on;
    plot(1:10,bfb0(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{2}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF B:B0 in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print('-dpng',[figp cfile3 '_allRE_nmort2.png'])

b0 = NaN*ones(length(frate),length(RE),length(spots));
for f=1:length(frate)
    b0(f,:,:) = B0(:,:,3);
end
bfb0 = m3_BF ./ b0;
for s = 1:length(spots)
    figure(8);
    subplot(4,3,s)
    plot(1:10,y,'--k'); hold on;
    plot(1:10,bfb0(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{3}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF B:B0 in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print('-dpng',[figp cfile3 '_allRE_nmort3.png'])

b0 = NaN*ones(length(frate),length(RE),length(spots));
for f=1:length(frate)
    b0(f,:,:) = B0(:,:,4);
end
bfb0 = m4_BF ./ b0;
for s = 1:length(spots)
    figure(9);
    subplot(4,3,s)
    plot(1:10,y,'--k'); hold on;
    plot(1:10,bfb0(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{4}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF B:B0 in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print('-dpng',[figp cfile3 '_allRE_nmort4.png'])

b0 = NaN*ones(length(frate),length(RE),length(spots));
for f=1:length(frate)
    b0(f,:,:) = B0(:,:,5);
end
bfb0 = m5_BF ./ b0;
for s = 1:length(spots)
    figure(10);
    subplot(4,3,s)
    plot(1:10,y,'--k'); hold on;
    plot(1:10,bfb0(:,:,s),'LineWidth',1.5)
    xlim([0 length(frate)+1])
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{5}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF B:B0 in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg,'LineWidth',1.5)
xlim([0 length(frate)+1])
ylim([0 0.1])
legend(sreff)
stamp(cfile2)
print('-dpng',[figp cfile3 '_allRE_nmort5.png'])

%% Plot Bmsy:B0
cfile4 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_BE' BE '_CC' CC '_MF_BmsyB0'];

[maxD0,fid] = max(m0_maxC,[],1);
fid = squeeze(fid);
bmsy = NaN*ones(length(RE),length(spots));
for r = 1:length(RE)
    for s = 1:length(spots)
        bmsy(r,s) = m0_BF(fid(r,s),r,s);
    end
end
bb = bmsy./B0(:,:,1);
for s = 1:length(spots)
    figure(11);
    subplot(4,3,s)
    plot(0:9,y,'--k'); hold on;
    plot(1:3,bb(:,s),'.k','MarkerSize',15)
    xlim([0 length(RE)+1])
    set(gca,'XTick',1:length(RE),'XTickLabel',reff);
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{1}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF Bmsy:B0 in final year')
    end
end
stamp(cfile2)
print('-dpng',[figp cfile4 '_allRE_nmort0.png'])

[maxD2,fid] = max(m2_maxC,[],1);
fid = squeeze(fid);
bmsy = NaN*ones(length(RE),length(spots));
for r = 1:length(RE)
    for s = 1:length(spots)
        bmsy(r,s) = m2_BF(fid(r,s),r,s);
    end
end
bb = bmsy./B0(:,:,2);
for s = 1:length(spots)
    figure(12);
    subplot(4,3,s)
    plot(0:9,y,'--k'); hold on;
    plot(1:3,bb(:,s),'.k','MarkerSize',15)
    xlim([0 length(RE)+1])
    set(gca,'XTick',1:length(RE),'XTickLabel',reff);
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{2}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF Bmsy:B0 in final year')
    end
end
stamp(cfile2)
print('-dpng',[figp cfile4 '_allRE_nmort2.png'])

[maxD3,fid] = max(m3_maxC,[],1);
fid = squeeze(fid);
bmsy = NaN*ones(length(RE),length(spots));
for r = 1:length(RE)
    for s = 1:length(spots)
        bmsy(r,s) = m3_BF(fid(r,s),r,s);
    end
end
bb = bmsy./B0(:,:,3);
for s = 1:length(spots)
    figure(13);
    subplot(4,3,s)
    plot(0:9,y,'--k'); hold on;
    plot(1:3,bb(:,s),'.k','MarkerSize',15)
    xlim([0 length(RE)+1])
    set(gca,'XTick',1:length(RE),'XTickLabel',reff);
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{3}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF Bmsy:B0 in final year')
    end
end
stamp(cfile2)
print('-dpng',[figp cfile4 '_allRE_nmort3.png'])

[maxD4,fid] = max(m4_maxC,[],1);
fid = squeeze(fid);
bmsy = NaN*ones(length(RE),length(spots));
for r = 1:length(RE)
    for s = 1:length(spots)
        bmsy(r,s) = m4_BF(fid(r,s),r,s);
    end
end
bb = bmsy./B0(:,:,4);
for s = 1:length(spots)
    figure(14);
    subplot(4,3,s)
    plot(0:9,y,'--k'); hold on;
    plot(1:3,bb(:,s),'.k','MarkerSize',15)
    xlim([0 length(RE)+1])
    set(gca,'XTick',1:length(RE),'XTickLabel',reff);
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{4}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF Bmsy:B0 in final year')
    end
end
stamp(cfile2)
print('-dpng',[figp cfile4 '_allRE_nmort4.png'])

[maxD5,fid] = max(m5_maxC,[],1);
fid = squeeze(fid);
bmsy = NaN*ones(length(RE),length(spots));
for r = 1:length(RE)
    for s = 1:length(spots)
        bmsy(r,s) = m5_BF(fid(r,s),r,s);
    end
end
bb = bmsy./B0(:,:,5);
for s = 1:length(spots)
    figure(15);
    subplot(4,3,s)
    plot(0:9,y,'--k'); hold on;
    plot(1:3,bb(:,s),'.k','MarkerSize',15)
    xlim([0 length(RE)+1])
    set(gca,'XTick',1:length(RE),'XTickLabel',reff);
    ylim([0 1.5])
    if (s==2)
        str = {['M=' Mort{5}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('MF Bmsy:B0 in final year')
    end
end
stamp(cfile2)
print('-dpng',[figp cfile4 '_allRE_nmort5.png'])


%% Find best
[Rmesh,Fmesh] = meshgrid(reff,Frate);
Cmax = nan*ones(length(spots),2,length(nmrt));
for s=1:11
    maxD0 = m0_maxC(:,:,s);
    m0 = find(maxD0==max(maxD0(:)));
    Cmax(s,1,1) = Rmesh(m0);
    Cmax(s,2,1) = Fmesh(m0);
    
    maxD2 = m2_maxC(:,:,s);
    m2 = find(maxD2==max(maxD2(:)));
    Cmax(s,1,2) = Rmesh(m2);
    Cmax(s,2,2) = Fmesh(m2);
    
    maxD3 = m3_maxC(:,:,s);
    m3 = find(maxD3==max(maxD3(:)));
    Cmax(s,1,3) = Rmesh(m3);
    Cmax(s,2,3) = Fmesh(m3);
    
    maxD4 = m4_maxC(:,:,s);
    m4 = find(maxD4==max(maxD4(:)));
    Cmax(s,1,4) = Rmesh(m4);
    Cmax(s,2,4) = Fmesh(m4);
    
    maxD5 = m5_maxC(:,:,s);
    m5 = find(maxD5==max(maxD5(:)));
    Cmax(s,1,5) = Rmesh(m5);
    Cmax(s,2,5) = Fmesh(m5);
end

%%
save([figp cfile2 '_allRE_allMort.mat'],'Cmax','m0_maxC','m2_maxC',...
    'm3_maxC','m4_maxC','m5_maxC','m0_BF','m2_BF','m3_BF','m4_BF','m5_BF',...
    'B0');

