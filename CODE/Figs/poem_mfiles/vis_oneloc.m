% Visualize output of POEM
% Spinup at one location
% 100 years

clear all
close all

dpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

sname = 'Spinup_NS_vspawn_';
loc = 'NS';
pi = csvread([dpath sname 'PISC.csv']);
pl = csvread([dpath sname 'PLAN.csv']);
de = csvread([dpath sname 'DETR.csv']);

%% Plots over time
x=1:length(pi);
x=x/365;

%% Piscivore
figure(1)
subplot(1,2,1)
plot(x,pi)
xlim([x(1) x(end)])
title(['Piscivore ' loc])
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
legend('1','2','3','4','5','6','7','8','9','10')
subplot(2,2,2)
plot(x(1:730),pi(1:730,:),'Linewidth',2)
xlim([x(1) x(730)])
title(['Piscivore ' loc])
subplot(2,2,4)
plot(x((end-731):end),pi((end-731):end,:),'Linewidth',2)
xlim([x(end-731) x(end)])
print('-dpng',[fpath sname 'oneloc_pisc_time.png'])

% Planktivore
figure(2)
subplot(1,2,1)
plot(x,pl)
xlim([x(1) x(end)])
title(['Planktivore ' loc])
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
legend('1','2','3','4','5','6','7','8','9','10')
subplot(2,2,2)
plot(x(1:730),pl(1:730,:),'Linewidth',2)
xlim([x(1) x(730)])
title(['Planktivore ' loc])
subplot(2,2,4)
plot(x((end-366):end),pl((end-366):end,:),'Linewidth',2)
xlim([x(end-366) x(end)])
print('-dpng',[fpath sname 'oneloc_plan_time.png'])

%Detritivore
figure(3)
subplot(1,2,1)
plot(x,de)
xlim([x(1) x(end)])
title(['Detritivore ' loc])
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
legend('1','2','3','4','5','6','7','8','9','10')
subplot(2,2,2)
plot(x(1:730),de(1:730,:),'Linewidth',2)
xlim([x(1) x(730)])
title(['Detritivore ' loc])
subplot(2,2,4)
plot(x((end-366):end),de((end-366):end,:),'Linewidth',2)
xlim([x(end-366) x(end)])
print('-dpng',[fpath sname 'oneloc_detr_time.png'])

% All size classes of all
figure(4)
for i=1:10
    subplot(3,4,i)
    plot(x(36134:36500),pi(36134:36500,i),'k','Linewidth',2); hold on
    xlim([x(36134) x(end)])
    %ylim([0 5e-3])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
end
%print('-dpng',[fpath sname 'oneloc_pisc_sizes.png'])
%
figure(5)
for i=1:10
    subplot(3,4,i)
    plot(x(36134:36500),pl(36134:36500,i),'b','Linewidth',2); hold on
    xlim([x(36134) x(end)])
    %ylim([0 2.5e-16])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
end
%print('-dpng',[fpath sname 'oneloc_plan_sizes.png'])
%
figure(6)
for i=1:10
    subplot(3,4,i)
    plot(x(36134:36500),de(36134:36500,i),'r','Linewidth',2); hold on
    xlim([x(36134) x(end)])
    %ylim([0 1.2e-8])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
end
%print('-dpng',[fpath sname 'oneloc_detr_sizes.png'])
%%
figure(7)
for i=1:9
    subplot(3,4,i)
    plot(x(36134:36500),log(pi(36134:36500,i)),'k','Linewidth',2); hold on
    plot(x(36134:36500),log(pl(36134:36500,i)),'b','Linewidth',2); hold on
    plot(x(36134:36500),log(de(36134:36500,i)),'r','Linewidth',2); hold on
    xlim([x(36134) x(end)])
    ylim([-100 0])
    xlabel('Time (y)')
    ylabel('log Biomass (g km^-^2)')
end
subplot('position',[0.372 0.12 0.305 0.21])
plot(x(36134:36500),log(pi(36134:36500,i)),'k','Linewidth',2); hold on
plot(x(36134:36500),log(pl(36134:36500,i)),'b','Linewidth',2); hold on
plot(x(36134:36500),log(de(36134:36500,i)),'r','Linewidth',2); hold on
xlim([x(36134) x(end)])
ylim([-100 0])
xlabel('Time (y)')
ylabel('log Biomass (g km^-^2)')
legend('Piscivore','Planktivore','Detritivore')
legend('location','eastoutside')
%print('-dpng',[fpath sname 'oneloc_all_sizes.png'])

%% Final mean biomass size spectrum
t=1:length(pi);
lyr=t((end-365):end);
pi_sum=sum(pi(lyr,:));
pi_mean=mean(pi(lyr,:));
pl_sum=sum(pl(lyr,:));
pl_mean=mean(pl(lyr,:));
de_sum=sum(de(lyr,:));
de_mean=mean(de(lyr,:));

figure(8)
subplot(2,3,1)
bar(pi_sum,'k')
xlim([0 11])
title('Piscivores')
ylabel('Total Annual Biomass (g km^-^2)')
subplot(2,3,4)
bar(pi_mean,'k')
xlim([0 11])
ylabel('Mean Annual Biomass (g km^-^2)')

subplot(2,3,2)
bar(pl_sum,'b')
xlim([0 11])
title({loc; 'Planktivores'})
xlabel('Size class')
subplot(2,3,5)
bar(pl_mean,'b')
xlim([0 11])
xlabel('Size class')

subplot(2,3,3)
bar(de_sum,'r')
xlim([0 11])
title('Detritivores')
subplot(2,3,6)
bar(de_mean,'r')
xlim([0 11])
print('-dpng',[fpath sname 'oneloc_all_biomass_spec.png'])

%% log scale with weight
%Number of size classes
PI_N=10;
PL_N=10;
DE_N=10;

%Min body size (g)
PI_smin = 10;
PL_smin = .1;
DE_smin = .1;

%Max body size (g)
PI_smax = 10000;
PL_smax = 500;
DE_smax = 500;

%Body mass linearly distributed (g)
PI_s = linspace((PI_smin),(PI_smax),PI_N);
PL_s = linspace((PL_smin),(PL_smax),PL_N);
DE_s = linspace((DE_smin),(DE_smax),DE_N);

figure(9)
subplot(2,1,1)
plot(log(PI_s),log(pi_sum),'k','Linewidth',2); hold on;
plot(log(PL_s),log(pl_sum),'b','Linewidth',2); hold on;
plot(log(DE_s),log(de_sum),'r','Linewidth',2); hold on;
xlabel('log Weight of size class (g)')
ylabel('log Total Annual Biomass (g km^-^2)')

subplot(2,1,2)
plot(log(PI_s),log(pi_mean),'k','Linewidth',2); hold on;
plot(log(PL_s),log(pl_mean),'b','Linewidth',2); hold on;
plot(log(DE_s),log(de_mean),'r','Linewidth',2); hold on;
xlabel('log Weight of size class (g)')
ylabel('log Mean Annual Biomass (g km^-^2)')
%print('-dpng',[fpath sname 'oneloc_all_logbiomass_spec.png'])

%
figure(10)
subplot(3,1,1)
plot(log(PI_s),log(pi_sum),'k','Linewidth',2); hold on;
xlim([-3 10])
title({loc; 'Piscivores'})

subplot(3,1,2)
plot(log(PL_s),log(pl_mean),'b','Linewidth',2); hold on;
xlim([-3 10])
ylabel('log Mean Annual Biomass (g km^-^2)')
title('Planktivores')

subplot(3,1,3)
plot(log(DE_s),log(de_mean),'r','Linewidth',2); hold on;
xlim([-3 10])
xlabel('log Weight of size class (g)')
title('Detritivores')
print('-dpng',[fpath sname 'oneloc_each_logbiomass_spec.png'])






