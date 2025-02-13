%Compare swimming speed of Megrey to a BL-based one

clear all
close all

load('/Volumes/GFDL/GCM_DATA/CORE-forced/ocean.186101-200512.temp_100_avg.mat',...
    'TEMP_100')
T100 = TEMP_100(:,:,end) - 273.15;
mid = find(T100 == min(T100(:)));
T100(mid) = NaN;

clear TEMP_100 nmdz_avg200_88 nlgz_avg200_88

%%
L_s = 10.0; % small
L_m = 200.0; % medium
L_l = 1.0e3;% large

%%! Mass from length using Andersen & Beyer 2013
% Convert from mm to cm and use their const coeff = 0.01g/cm3
M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;
    
wgt = M_m; 	
wL = M_l;

T = 15.0;
w = ((3.9*wgt.^0.13 * exp(0.149*T100)) /100);
Q = zeros(360,200);
aa = find(~isnan(T100));
Q(aa) = w(aa);

w2 = exp(0.063*(T100-15.0)) * 0.5*L_m*1e-3;
Q2 = zeros(360,200);
Q2(aa) = w2(aa);

w3 = ((3.9*M_l.^0.13 * exp(0.149*T100)) /100);
Q3 = zeros(360,200);
Q3(aa) = w3(aa);

w4 = exp(0.063*(T100-15.0)) * 0.5*L_l*1e-3;
Q4 = zeros(360,200);
Q4(aa) = w4(aa);

%%
figure(1)
pcolor(w')
shading flat
colorbar
title('Megrey medium')

figure(2)
pcolor(w2')
shading flat
colorbar
title('half BL medium')

figure(3)
pcolor(w3')
shading flat
colorbar
title('Megrey large')

figure(4)
pcolor(w4')
shading flat
colorbar
title('half BL large')

%%
%t=unique(T100);
t = min(T100(:)):max(T100(:));
swim1=((3.9*wgt.^0.13 * exp(0.149*t)) /100);
swim2 = exp(0.063*(t-15.0)) * 0.5*L_m*1e-3;
swim3=((3.9*wL.^0.13 * exp(0.149*t)) /100);
swim4 = exp(0.063*(t-15.0)) * 0.5*L_l*1e-3;

%%
figure(5)
subplot(2,1,1)
plot(t,swim1,'b','LineWidth',2); hold on;
plot(t,swim2,'r','LineWidth',2);
legend('Megrey-wgt','halfBL')
legend('location','northeast')
ylim([0 2*L_m*1e-3])
xlim([-2 32])
title('Medium')
ylabel('speed (m/s)')

subplot(2,1,2)
plot(t,swim3,'b','LineWidth',2); hold on;
plot(t,swim4,'r','LineWidth',2);
ylim([0 2*L_l*1e-3])
xlim([-2 32])
title('Large')
ylabel('speed (m/s)')
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_other/biol_eqs/swim_speed_alt_big_sizes.png')
