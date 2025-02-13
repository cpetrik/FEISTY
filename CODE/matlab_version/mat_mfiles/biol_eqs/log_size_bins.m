% Test log spaced size bins
% Then calc Maturity

clear all
close all

%% Current
PI_smin = 10;   %min size in g
PI_smax = 1e4;  %max size in g
PI_N = 10;      %number of size classes
lPI_s = linspace(log10(PI_smin),log10(PI_smax),PI_N);
PI_s = 10.^(lPI_s);
ds_pi = abs(PI_s - (PI_s(end) * 0.1));
pid = find(ds_pi == min(ds_pi));
pid = pid(1);
PI_K = [ones(1,pid), linspace(1,0,PI_N - pid)];
PI_rep = 1-PI_K;

figure(1)
subplot(1,2,1)
plot(PI_K,'LineWidth',2);
axis([1 10 0 1.1])
xlabel('Size class')
ylabel('Kappa = fraction growth')
title('Piscivore')
subplot(1,2,2)
plot(PI_rep,'LineWidth',2);
axis([1 10 0 1.1])
xlabel('Size class')
ylabel('1-Kappa = fraction repro')
title('Piscivore')

%% Andersen
nm = 0.25;
n=0.75;
psi = (1+(PI_s/(nm.*PI_s(end))).^-10).^-1 .* (PI_s/PI_s(end)).^(1-n);
grow=1-psi;

psiA = psi;
growA = grow;

figure(2)
subplot(1,2,1)
plot(growA,'LineWidth',2);
axis([1 10 0 1.1])
xlabel('Size class')
ylabel('1-psi = immature -> growth')
title('Piscivore')
subplot(1,2,2)
plot(psiA,'LineWidth',2);
axis([1 10 0 1.1])
xlabel('Size class')
ylabel('psi = mature -> repro')
title('Piscivore')

%% My FB & RAM analysis of Lm vs Linf
nm = 0.1213;
n=0.75;
psi = (1+(PI_s/(nm.*PI_s(end))).^-10).^-1 .* (PI_s/PI_s(end)).^(1-n);
grow=1-psi;

figure(3)
subplot(1,2,1)
plot(grow,'LineWidth',2);
axis([1 10 0 1.1])
xlabel('Size class')
ylabel('1-psi = immature -> growth')
title('Piscivore')
subplot(1,2,2)
plot(psi,'LineWidth',2);
axis([1 10 0 1.1])
xlabel('Size class')
ylabel('psi = mature -> repro')
title('Piscivore')

%% All together
figure(4)
subplot(1,2,1)
plot(PI_s,PI_K,'k','LineWidth',2); hold on;
plot(PI_s,growA,'b','LineWidth',2); hold on;
plot(PI_s,grow,'r','LineWidth',2); hold on;
ylim([0 1.1])
legend('POEM','Andersen eta','My eta')
xlabel('Weight (g)')
ylabel('Fraction growth/immature')
title('Piscivore')
subplot(1,2,2)
plot(PI_s,PI_rep,'k','LineWidth',2);hold on;
plot(PI_s,psiA,'b','LineWidth',2); hold on;
plot(PI_s,psi,'r','LineWidth',2); hold on;
ylim([0 1.1])
xlabel('Weight (g)')
ylabel('Fraction repro/mature')
title('Piscivore')
print('-dpng','maturity_fns_pisc_logsize.png')


