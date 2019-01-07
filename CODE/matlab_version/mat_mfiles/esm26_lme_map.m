figure(10)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,tlme)
colormap('prism')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('ESM2.6 1^o LMEs')
print('-dpng',[cpath 'esm26_LME.png'])

figure(100)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,tlme)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 66]);
hcb = colorbar('h');
ylim(hcb,[0 66])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('ESM2.6 1^o LMEs')
print('-dpng',[cpath 'esm26_LME_num.png'])
