%Compare swimming speed of Megrey to a BL-based one

clear 
close all

%%
load('/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/ocean.186101-200512.temp_100_avg.mat',...
    'TEMP_100')
T100 = TEMP_100(:,:,end) - 273.15;
mid = find(T100 == min(T100(:)));
T100(mid) = NaN;

clear TEMP_100 

%%
%!Individual Mass (g) = geometric mean
M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small   4.6-36.8 mm, 13.1 mm
L_m = 10.0 * (M_m/0.01)^(1/3); % medium 36.8-292 mm, 10.4 cm
L_l = 10.0 * (M_l/0.01)^(1/3); % large  0.292-2.32 m, 0.82 m
    
%%
wgt = M_m; 	
wL = M_l;

T = 15.0;
%w = ((3.9*wgt.^0.13 * exp(0.149*T100)) /100); %Megrey
w = ((3.9*wgt.^0.13 * exp(0.063*(T100-15.0))) /100); %Megrey wgt-dep, COBALT T-dep
Q = zeros(360,200);
aa = (~isnan(T100));
Q(aa) = w(aa);

w2 = exp(0.063*(T100-15.0)) * 0.5*L_m*1e-3;
Q2 = zeros(360,200);
Q2(aa) = w2(aa);

%w3 = ((3.9*M_l.^0.13 * exp(0.149*T100)) /100); %Megrey
w3 = ((3.9*M_l.^0.13 * exp(0.063*(T100-15.0))) /100); %Megrey wgt-dep, COBALT T-dep
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
%swim1=((3.9*wgt.^0.13 * exp(0.149*t)) /100);
swim1=((3.9*wgt.^0.13 * exp(0.063*(t-15.0))) /100);
swim2 = exp(0.063*(t-15.0)) * 0.5*L_m*1e-3;
%swim3=((3.9*wL.^0.13 * exp(0.149*t)) /100);
swim3=((3.9*wL.^0.13 * exp(0.063*(t-15.0))) /100);
swim4 = exp(0.063*(t-15.0)) * 0.5*L_l*1e-3;

%%
figure(5)
subplot(2,1,1)
plot(t,swim1,'b','LineWidth',2); hold on;
plot(t,swim2,'r','LineWidth',2);
legend('Megrey-wgt','halfBL')
legend('location','northwest')
%ylim([0 2*L_m*1e-3])
xlim([-2 32])
title('Medium')
ylabel('speed (m/s)')

subplot(2,1,2)
plot(t,swim3,'b','LineWidth',2); hold on;
plot(t,swim4,'r','LineWidth',2);
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
title('Large')
ylabel('speed (m/s)')
%print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_other/biol_eqs/swim_speed_alt_big_sizes.png')

%% Megrey wgt-dep, COBALT T-dep
blM = swim1 / (L_m*1e-3);
blL = swim3 / (L_l*1e-3);

figure(6)
subplot(2,1,1)
plot(t,swim1,'r','LineWidth',2); hold on;
plot(t,swim3,'b','LineWidth',2); hold on;
legend('M','L')
legend('location','northwest')
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
title('Megrey wgt-dep, COBALT T-dep')
ylabel('speed (m/s)')

subplot(2,1,2)
plot(t,blM,'r','LineWidth',2); hold on;
plot(t,blL,'b','LineWidth',2); hold on;
legend('M','L')
legend('location','northwest')
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
ylabel('speed (BL/s)')

%% Megrey wgt-dep, COBALT T-dep, but decr wgt-dep constant so uM < 1 BL/s
swimM=((2.02*M_m.^0.13 * exp(0.063*(t-10.0))) /100);
swimL=((2.02*M_l.^0.13 * exp(0.063*(t-10.0))) /100);
blM = swimM / (L_m*1e-3);
blL = swimL / (L_l*1e-3);

figure(7)
subplot(2,1,1)
plot(t,swimM,'r','LineWidth',2); hold on;
plot(t,swimL,'b','LineWidth',2); hold on;
legend('M','L')
legend('location','northwest')
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
title('Megrey wgt-dep, COBALT T-dep')
ylabel('speed (m/s)')

subplot(2,1,2)
plot(t,blM,'r','LineWidth',2); hold on;
plot(t,blL,'b','LineWidth',2); hold on;
legend('M','L')
legend('location','northwest')
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
ylabel('speed (BL/s)')
