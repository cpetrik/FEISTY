y=time;

figure
subplot(3,3,1)
plot(y,(sf_tmean),'Linewidth',1); hold on;
ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,2)
plot(y,(mf_tmean),'Linewidth',1); hold on;
ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,3)
plot(y,(b_tmean),'Linewidth',1); hold on;
%ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,4)
plot(y,(sp_tmean),'Linewidth',1); hold on;
ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,5)
plot(y,(mp_tmean),'Linewidth',1); hold on;
ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,6)
plot(y,(lp_tmean),'Linewidth',1); hold on;
ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,7)
plot(y,(sd_tmean),'Linewidth',1); hold on;
ylim([0 1e-6])
xlim([1850 1860])

subplot(3,3,8)
plot(y,(md_tmean),'Linewidth',1); hold on;
ylim([0 1e-5])
xlim([1850 1860])

subplot(3,3,9)
plot(y,(ld_tmean),'Linewidth',1); hold on;
%ylim([0 1e-6])
xlim([1850 1860])

%%
figure
subplot(3,1,1)
plot(y,(b_tmean),'Linewidth',1); hold on;
%ylim([0 1e-5])
%xlim([1850 1860])

subplot(3,1,2)
plot(y,(md_tmean),'Linewidth',1); hold on;
%ylim([0 1e-5])
%xlim([1850 1860])

subplot(3,1,3)
plot(y,(ld_tmean),'Linewidth',1); hold on;
%ylim([0 1e-6])
%xlim([1850 1860])