% Visualize advection test cases w/diffusion
% Test 1st order upwind, 2nd order centered, 4th order centered
% advection schemes from MOM5
% MOM6 grid

clear all
close all

dpath = '/Volumes/GFDL/CSV/advect_tests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/Advec-Diff/figs_tests/';

cname = 'AtlArc_even_dt100hr_mom6_preind_vel_surf';
load([dpath 'Matlab_adv_diff_1st_2nd_4th_K0_' cname '_UV2.mat'],'B1','B2','B4');

load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid1D_MOM6_preindust.mat','GRD');
ID = GRD.ID;
clear GRD
load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid2D_MOM6_preindust.mat')

%%
bio1 = B1;
bio2 = B2;
bio4 = B4;

clear B1 B2 B4

%% Conservation of mass
% 1st order upwind
[ni,nj,nd] = size(bio1);

mass1 = bio1 .* repmat(GRD.AREA,1,1,nd);
totb1 = squeeze(nansum(nansum(mass1,1)));
cons1 = 100*(totb1(end)-totb1(1))/totb1(1)

%
yrs=[1:length(totb1)]/365;
figure(51)
plot(yrs,totb1,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles 1st')
print('-dpng',[fpath 'advec_diff_test_K0_' cname '_UV2_totb_1st.png'])

%% 2nd order
mass2 = bio2 .* repmat(GRD.AREA,1,1,nd);
totb2 = squeeze(nansum(nansum(mass2,1)));
cons2 = 100*(totb2(end)-totb2(1))/totb2(1)

%
yrs=[1:length(totb2)]/365;
figure(52)
plot(yrs,totb2,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles 2nd')
print('-dpng',[fpath 'advec_diff_test_' cname '_UV2_totb_2nd.png'])

%% 4th order
mass4 = bio4 .* repmat(GRD.AREA,1,1,nd);
totb4 = squeeze(nansum(nansum(mass4,1)));
cons4 = 100*(totb4(end)-totb4(1))/totb4(1)

%
yrs=[1:length(totb4)]/365;
figure(54)
plot(yrs,totb4,'LineWidth',2)
%xlim([0 12])
%xlabel('days')
xlabel('Year')
%title(cname)
ylabel('Total number of particles 4th')
print('-dpng',[fpath 'advec_diff_test_' cname '_UV2_totb_4th.png'])

%% plot info
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

% Land
surf_tmask = GRD.lmask;
lmask = surf_tmask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

%axes labels
xt = -250:50:50;
xl = xt;
xl(xl<-180) = xl(xl<-180) + 350;

t = 1:72.75:nd;
%t = 1:3:15;
t = round(t);

% colors
cmB=cbrewer('seq','Blues',50);

%% Global flat
for n=1:length(t)
    B1 = bio1(:,:,t(n));
%     B2 = bio2(:,:,t(n));
%     B4 = bio4(:,:,t(n));
    
    figure(1)
    clf
    axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
        'Grid','off','FLineWidth',1,'origin',[0 -100 0])
    surfm(GRD.LAT,GRD.LON,B1)
    load coast;                     %decent looking coastlines
    h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
    colorbar;
    caxis([0 200]);
    colormap(cmB)
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_1st_K0_' cname '_UV2_' num2str(t(n)) '.png'])
    
%     figure(2)
%     clf
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%     surfm(GRD.LAT,GRD.LON,B2)
%     load coast;                     %decent looking coastlines
%     h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     colorbar;
%     caxis([0 200]);
%     colormap(cmB)
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[fpath 'advec_diff_test_2nd_' cname '_UV2_' num2str(t(n)) '.png'])
%     
%     figure(4)
%     clf
%     axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%         'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%     surfm(GRD.LAT,GRD.LON,B4)
%     load coast;                     %decent looking coastlines
%     h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     colorbar;
%     caxis([0 200]);
%     colormap(cmB)
%     title(['Day ' num2str(t(n)) ' Year 1'])
%     print('-dpng',[fpath 'advec_diff_test_4th_' cname '_UV2_' num2str(t(n)) '.png'])
end


%% Arctic projection
for n=1:length(t)
    B1 = bio1(:,:,t(n));
    %     B2 = bio2(:,:,t(n));
    %     B4 = bio4(:,:,t(n));
    
    figure(1)
    clf
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(GRD.LON,GRD.LAT,B1);
    shading interp
    colorbar
    %colormap('jet')
    colormap(cmB)
    caxis([0 110])
    m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_1st_K0_' cname '_UV2_arcticproj_' num2str(t(n)) '.png'])
    
    %     figure(2)
    %     clf
    %     m_proj('stereographic','lat',90,'long',30,'radius',30);
    %     m_pcolor(GRD.LON,GRD.LAT,B2);
    %     shading interp
    %     colorbar
    %     %colormap('jet')
    %     colormap(cmB)
    %     caxis([0 150])
    %     m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    %     m_coast('patch',[.7 .7 .7],'edgecolor','k');
    %     title(['Day ' num2str(t(n)) ' Year 1'])
    %     print('-dpng',[fpath 'advec_diff_test_2nd_K100_' cname '_arcticproj_' num2str(t(n)) '.png'])
    %
    %     figure(4)
    %     clf
    %     m_proj('stereographic','lat',90,'long',30,'radius',30);
    %     m_pcolor(GRD.LON,GRD.LAT,B4);
    %     shading interp
    %     colorbar
    %     %colormap('jet')
    %     colormap(cmB)
    %     caxis([0 150])
    %     m_grid('xtick',6,'tickdir','out','ytick',[70 80],'linest','-');
    %     m_coast('patch',[.7 .7 .7],'edgecolor','k');
    %     title(['Day ' num2str(t(n)) ' Year 1'])
    %     print('-dpng',[fpath 'advec_diff_test_4th_K100_' cname '_arcticproj_' num2str(t(n)) '.png'])
    
end

%% Antarctic projection
for n=1:length(t)
    B1 = bio1(:,:,t(n));
    
    figure(2)
    clf
    m_proj('stereographic','lat',-90,'long',30,'radius',50);
    m_pcolor(GRD.LON,GRD.LAT,B1);
    shading interp
    colorbar
    %colormap('jet')
    colormap(cmB)
    caxis([0 100])
    m_grid('xtick',12,'tickdir','out','ytick',[-50 -60 -70],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(t(n)) ' Year 1'])
    print('-dpng',[fpath 'advec_diff_test_' cname '_Spoleproj_' num2str(t(n)) '.png'])
end

%% Movie

% figure(1)
% clf
% fig1=figure(1);
% %
% winsize = get(fig1,'Position');
% winsize(1:2) = [0 0];
% numframes=365;
% A=moviein(numframes,fig1,winsize);
% set(fig1,'NextPlot','replacechildren')
%
% % draw background
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1,'origin',[0 -100 0])
% load coast;                     %decent looking coastlines
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%
% %
% j=0;
% for d=1:365
%     j=j+1
%     surfm(GRD.LAT,GRD.LON,bio1(:,:,d))
%     colorbar;
%     caxis([0 150]);
%     colormap(cmB)
%     h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%     title(['Day ' num2str(d) ' Year 1'])
%     A(:,j)=getframe(fig1,winsize);
% end
%
% %
% %movie(fig1,A,1,3,winsize);
% save([dpath 'advec_test_1st_' cname '_UV2_movie.mat'], 'A')
%
% v1 = VideoWriter([fpath 'advec_test_1st_' cname '_UV2_movie.avi']);
% open(v1)
% writeVideo(v1,A);
% close(v1)
%
% % v2 = VideoWriter([fpath 'advec_test_1st_' cname '_UV2_movie_nocomp.avi'],...
% %     'Uncompressed AVI');
% % open(v2)
% % writeVideo(v2,A);
% % close(v2)

