i=1;
alt1 = alt{i};
    spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
    fpath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
        'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
    load([fpath 'grid_mortality_all_',alt1,'.mat'])
AfmortD = fmortD;
AfmortF = fmortF;
AfmortP = fmortP;
clear fmortD fmortF fmortP

i=2;
alt1 = alt{i};
    spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
    fpath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
        'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
    load([fpath 'grid_mortality_all_',alt1,'.mat'])
EfmortD = fmortD;
EfmortF = fmortF;
EfmortP = fmortP;
clear fmortD fmortF fmortP

i=3;
alt1 = alt{i};
    spath = ['/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/',alt1,'/'];
    fpath = ['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/',...
        'FEISTY_other/fishing_ms_ideas/fishing_effort_ms/fishing_for_FEISTY/',alt1,'/'];
    load([fpath 'grid_mortality_all_',alt1,'.mat'])
NfmortD = fmortD;
NfmortF = fmortF;
NfmortP = fmortP;

clear fmortD fmortF fmortP
%%
figure
subplot(3,3,1)
pcolor(AfmortD(28600:28650,1:50)); shading flat;
title('AD')
colorbar
caxis([0 0.05])

subplot(3,3,2)
pcolor(AfmortF(28600:28650,1:50)); shading flat;
title('AF')
colorbar
caxis([0 0.05])

subplot(3,3,3)
pcolor(AfmortP(28600:28650,1:50)); shading flat;
title('AP')
colorbar
caxis([0 0.2])

subplot(3,3,4)
pcolor(EfmortD(28600:28650,1:50)); shading flat;
title('ED')
colorbar
caxis([0 0.05])

subplot(3,3,5)
pcolor(EfmortF(28600:28650,1:50)); shading flat;
title('EF')
colorbar
caxis([0 0.05])

subplot(3,3,6)
pcolor(EfmortP(28600:28650,1:50)); shading flat;
title('EP')
colorbar
caxis([0 0.2])

subplot(3,3,7)
pcolor(NfmortD(28600:28650,1:50)); shading flat;
title('ND')
colorbar
caxis([0 0.05])

subplot(3,3,8)
pcolor(NfmortF(28600:28650,1:50)); shading flat;
title('NF')
colorbar
caxis([0 0.05])

subplot(3,3,9)
pcolor(NfmortP(28600:28650,1:50)); shading flat;
title('NP')
colorbar
caxis([0 0.2])

%%
testD=(AfmortD==NfmortD);
testD2= (AfmortD~=NfmortD);
sum(testD(:))
sum(testD2(:))

testF=(AfmortF==NfmortF);
testF2= (AfmortF~=NfmortF);
sum(testF(:))
sum(testF2(:))

testP=(AfmortP==NfmortP);
testP2= (AfmortP~=NfmortP);
sum(testP(:))
sum(testP2(:))

