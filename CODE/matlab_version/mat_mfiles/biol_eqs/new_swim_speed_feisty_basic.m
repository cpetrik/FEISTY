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
% temp-dep like COBALT
% wgt-dep like cmax or enc?
% power on enc fn 0.20
% power on cmax fn 0.25

Tref = 10.0;
aa = (~isnan(T100));
t = min(T100(:)):max(T100(:));

%Megrey wgt-dep, COBALT T-dep
wM = ((3.9*M_m.^0.13 * exp(0.063*(T100-Tref))) /100); 
QM = zeros(360,200);
QM(aa) = wM(aa);

wL = ((3.9*M_l.^0.13 * exp(0.063*(T100-Tref))) /100); 
QL = zeros(360,200);
QL(aa) = wL(aa);


%Megrey wgt-dep, COBALT T-dep, but decr wgt-dep constant so uM < 1 BL/s
wM2=((2.02*M_m.^0.13 * exp(0.063*(T100-Tref))) /100);
wL2=((2.02*M_l.^0.13 * exp(0.063*(T100-Tref))) /100);
QM2 = zeros(360,200);
QM2(aa) = wM2(aa);
QL2 = zeros(360,200);
QL2(aa) = wL2(aa);


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
% Megrey wgt-dep, COBALT T-dep, but decr wgt-dep constant so uM < 1 BL/s
swimM=((2.02*M_m.^0.13 * exp(0.063*(t-10.0))) /100);
swimL=((2.02*M_l.^0.13 * exp(0.063*(t-10.0))) /100);
blM = swimM / (L_m*1e-3);
blL = swimL / (L_l*1e-3);

% Cmax wgt-dep, COBALT T-dep
swimM2=((2.02*M_m.^0.25 * exp(0.063*(t-10.0))) /100);
swimL2=((2.02*M_l.^0.25 * exp(0.063*(t-10.0))) /100);
blM2 = swimM2 / (L_m*1e-3);
blL2 = swimL2 / (L_l*1e-3);

% Enc wgt-dep, COBALT T-dep
swimM3=((2.02*M_m.^0.2 * exp(0.063*(t-10.0))) /100);
swimL3=((2.02*M_l.^0.2 * exp(0.063*(t-10.0))) /100);
blM3 = swimM3 / (L_m*1e-3);
blL3 = swimL3 / (L_l*1e-3);


%%
figure(7)
subplot(2,1,1)
plot(t,swimM,'r','LineWidth',2); hold on;
plot(t,swimL,'b','LineWidth',2); hold on;
plot(t,swimM2,'m','LineWidth',2); hold on;
plot(t,swimL2,'c','LineWidth',2); hold on;
plot(t,swimM3,'color',[0.75 0 0.5],'LineWidth',2); hold on;
plot(t,swimL3,'color',[0 0.5 0.75],'LineWidth',2); hold on;
legend('M','L','M2','L2','M3','L3')
legend('location','northwest')
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
title('Megrey wgt-dep, COBALT T-dep')
ylabel('speed (m/s)')

subplot(2,1,2)
plot(t,blM,'r','LineWidth',2); hold on;
plot(t,blL,'b','LineWidth',2); hold on;
plot(t,blM2,'m','LineWidth',2); hold on;
plot(t,blL2,'c','LineWidth',2); hold on;
plot(t,blM3,'color',[0.75 0 0.5],'LineWidth',2); hold on;
plot(t,blL3,'color',[0 0.5 0.75],'LineWidth',2); hold on;
legend('M','L','M2','L2','M3','L3')
legend('location','northwest')
%ylim([0 2*L_l*1e-3])
xlim([-2 32])
ylabel('speed (BL/s)')
