% Visualize output of POEM
% Spinup at one location
% 50 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

sname = 'Spinup_GB_';
sname2 = 'Spinup_';
loc = 'GB';
SP = csvread([dpath sname 'Sml_p.csv']);
SF = csvread([dpath sname 'Sml_f.csv']);
SD = csvread([dpath sname 'Sml_d.csv']);
MP = csvread([dpath sname 'Med_p.csv']);
MF = csvread([dpath sname 'Med_f.csv']);
MD = csvread([dpath sname 'Med_d.csv']);
LP = csvread([dpath sname 'Lrg_p.csv']);

rMF = csvread([dpath sname2 'Rep_' loc '_Med_f.csv']);
rLP = csvread([dpath sname2 'Rep_' loc '_Lrg_p.csv']);
rMF = rMF .* MF;
rLP = rLP .* LP;

mSP = csvread([dpath sname2 'Rec_' loc '_Sml_p.csv']);
mSF = csvread([dpath sname2 'Rec_' loc '_Sml_f.csv']);
mMP = csvread([dpath sname2 'Rec_' loc '_Med_p.csv']);
mMF = csvread([dpath sname2 'Rec_' loc '_Med_f.csv']);
mLP = csvread([dpath sname2 'Rec_' loc '_Lrg_p.csv']);

%% Plots over time
x=1:length(SP);


%% Piscivore
figure(1)
plot(x,SP,'b','Linewidth',2); hold on;
plot(x,MP,'r','Linewidth',2); hold on;
plot(x,LP,'k','Linewidth',2); hold on;
xlim([x(1) x(end)])
title(['Pisc ' loc])
legend('S','M','L')

% Planktivore
figure(2)
plot(x,SF,'b','Linewidth',2); hold on;
plot(x,MF,'r','Linewidth',2); hold on;
xlim([x(1) x(end)])
title(['Forage fishes ' loc])
legend('S','M')

% Detritivore
figure(3)
plot(x,SD,'b','Linewidth',2); hold on;
plot(x,MD,'r','Linewidth',2); hold on;
xlim([x(1) x(end)])
title(['Demersal ' loc])
legend('S','M')


%% BIOMASS, REPRODUCTION, MATURATION

% Piscivore
figure(10)
subplot(3,1,1)
plot(x,SP,'b','Linewidth',2); hold on;
plot(x,MP,'r','Linewidth',2); hold on;
plot(x,LP,'k','Linewidth',2); hold on;
title(['Pelagic Piscivore ' loc])
legend('Larvae','Juveniles','Adults')

subplot(3,1,2)
plot(x,rLP,'k','Linewidth',2); hold on;
title('Repro')

subplot(3,1,3)
plot(x,mSP,'b','Linewidth',2); hold on;
plot(x,mMP,'r','Linewidth',2); hold on;
plot(x,mLP,'k','Linewidth',2); hold on;
title('Matur')

%%
figure(11)
subplot(3,1,1)
plot(x,SP,'b','Linewidth',2); hold on;
plot(x,mSP,'--m','Linewidth',2); hold on;
plot(x,rLP,':k','Linewidth',2); hold on;
legend('Larvae','A->L','produced')
xlim([1 10])

subplot(3,1,2)
plot(x,SP,'b','Linewidth',2); hold on;
plot(x,mMP,'--c','Linewidth',2); hold on;
plot(x,MP,'r','Linewidth',2); hold on;
legend('Larvae','L->J','Juveniles')
xlim([1 10])

subplot(3,1,3)
plot(x,MP,'r','Linewidth',2); hold on;
plot(x,mLP,'--','color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(x,LP,'k','Linewidth',2); hold on;
legend('Juveniles','J->A','Adults')
xlim([1 10])

%%
dsp=0;
dmp=0;
dlp=0;
dsp(2:365) = diff(SP);
dmp(2:365) = diff(MP);
dlp(2:365) = diff(LP);

figure(12)
subplot(3,1,1)
plot(x,dsp,'b','Linewidth',2); hold on;
plot(x,mSP,'--m','Linewidth',2); hold on;
plot(x,rLP,':k','Linewidth',2); hold on;
legend('Larvae','A->L','produced')
xlim([1 10])

subplot(3,1,2)
plot(x,dsp,'b','Linewidth',2); hold on;
plot(x,mMP,'--c','Linewidth',2); hold on;
plot(x,dmp,'r','Linewidth',2); hold on;
legend('Larvae','L->J','Juveniles')
xlim([1 10])

subplot(3,1,3)
plot(x,dmp,'r','Linewidth',2); hold on;
plot(x,mLP,'--','color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(x,dlp,'k','Linewidth',2); hold on;
legend('Juveniles','J->A','Adults')
xlim([1 10])


%% Planktivore
