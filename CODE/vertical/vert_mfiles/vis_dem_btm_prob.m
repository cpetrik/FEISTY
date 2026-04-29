
ld1 = squeeze(Ld_B(:,:,20:35,1));
figure(1)
for k=1:16
    subplot(4,4,k)
    pcolor(squeeze(ld1(:,:,k))); shading flat; colorbar;
end

%%
figure(2)
pcolor(squeeze(log10(ld1(:,:,14)+eps))); shading flat; colorbar;
%clim([-2 2])

figure(3)
pcolor(squeeze(log10(ld1(:,:,15)+eps))); shading flat; colorbar;

figure(4)
pcolor(squeeze(log10(ld1(:,:,16)+eps))); shading flat; colorbar;