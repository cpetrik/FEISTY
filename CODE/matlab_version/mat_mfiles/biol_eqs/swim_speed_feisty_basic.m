% Swimming speed
% Max recorded speeds are 56-110 km/h
% km/h -> m/s: *1e3 / 3600; -> 15.6-30.6 m/s

clear
close all

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

%% Watson et al. 2015
%u = kappa .* a .*s.^b; % s is mass in g
%a = 864;   % 864 m/d/g^0.13
%b = 0.13;
%kappa = 0.25-1.0 from temp-dep?
%m/d -> m/s
Uj_s = (864 *M_s.^0.13) /(60*60*24); % 864  /(60*60*24) = 0.01m/s (1 cm/s)
Uj_m = (864 *M_m.^0.13) /(60*60*24);
Uj_l = (864 *M_l.^0.13) /(60*60*24);

%% James' current 
%cm/s -> m/s
Uw_s = (3.9 *M_s.^0.13) * 1e-2; %/100*60*60*24;
Uw_m = (3.9 *M_m.^0.13) * 1e-2; %/100*60*60*24;
Uw_l = (3.9 *M_l.^0.13) * 1e-2; %/100*60*60*24;

%% Megrey et al with temp-dep (m/d) -> m/s
T = -2:0.2:37;
T1= -2:0.2:8.9;
T2= 9:0.2:37;
Ut_s = (1.0*M_s.^0.0 * exp(0.0*T)) /(60*60*24);
Ut_m1 = (3.9*M_m.^0.13 * exp(0.149*T1)) /(60*60*24);
Ut_m2 = (15.0*M_m.^0.13 * exp(0.0*T2)) /(60*60*24);
Ut_l1 = (3.9*M_l.^0.13 * exp(0.149*T1)) /(60*60*24);
Ut_l2 = (15.0*M_l.^0.13 * exp(0.0*T2)) /(60*60*24);
Ut_m=[Ut_m1 Ut_m2];
Ut_l=[Ut_l1 Ut_l2];

%% Q10 & BL (m/s)
%Tref = 34.5490; %gives 0.1 BL/s at -2C, but never above 1.3BL/s
%Tref = 22.5617; %gives 3 BL/s at 40C, but <1 BL/s below 22.5
Tref=18.5; %median temp over historical period
Ub_s = L_s * exp(0.063*(T-Tref)) *1e-3;
Ub_m = L_m * exp(0.063*(T-Tref)) *1e-3;
Ub_l = L_l * exp(0.063*(T-Tref)) *1e-3;

%% speed based on mass in m/s
figure(1)
plot(T,Ut_s,'b','LineWidth',2); hold on;
plot(T,Ut_m,'r','LineWidth',2); hold on;
plot(T,Ut_l,'k','LineWidth',2); hold on;
plot(T,Uw_s*ones(length(T),1),'c','LineWidth',2); hold on;
plot(T,Uw_m*ones(length(T),1),'m','LineWidth',2); hold on;
plot(T,Uw_l*ones(length(T),1),'color',[0.66 0.66 0.66],'LineWidth',2); hold on;
plot(T,Uj_s*ones(length(T),1),'color',[0 0.5 0.75],'LineWidth',2); hold on;
plot(T,Uj_m*ones(length(T),1),'color',[0.75 0 0.5],'LineWidth',2); hold on;
plot(T,Uj_l*ones(length(T),1),'color',[0.33 0.33 0.33],'LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL','JS','JM','JL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (m/s)')

%% speed based on mass in BL/s
figure(2)
plot(T,Ut_s/(1e-3*L_s),'b','LineWidth',2); hold on;
plot(T,Ut_m/(1e-3*L_m),'r','LineWidth',2); hold on;
plot(T,Ut_l/(1e-3*L_l),'k','LineWidth',2); hold on;
plot(T,Uw_s*ones(length(T),1)/(1e-3*L_s),'c','LineWidth',2); hold on;
plot(T,Uw_m*ones(length(T),1)/(1e-3*L_m),'m','LineWidth',2); hold on;
plot(T,Uw_l*ones(length(T),1)/(1e-3*L_l),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(T,Uj_s*ones(length(T),1)/(1e-3*L_s),'color',[0 0.5 0.75],'LineWidth',2); hold on;
plot(T,Uj_m*ones(length(T),1)/(1e-3*L_m),'color',[0.75 0 0.5],'LineWidth',2); hold on;
plot(T,Uj_l*ones(length(T),1)/(1e-3*L_l),'color',[0.33 0.33 0.33],'LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL','JS','JM','JL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (BL/s)')

%% speed based on mass & BL-Q10 in m/s
figure(3)
plot(T,Ut_s,'b','LineWidth',2); hold on;
plot(T,Ut_m,'r','LineWidth',2); hold on;
plot(T,Ut_l,'k','LineWidth',2); hold on;
plot(T,Uw_s*ones(length(T),1),'c','LineWidth',2); hold on;
plot(T,Uw_m*ones(length(T),1),'m','LineWidth',2); hold on;
plot(T,Uw_l*ones(length(T),1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(T,Uj_s*ones(length(T),1),'color',[0 0.5 0.75],'LineWidth',2); hold on;
plot(T,Uj_m*ones(length(T),1),'color',[0.75 0 0.5],'LineWidth',2); hold on;
plot(T,Uj_l*ones(length(T),1),'color',[0.33 0.33 0.33],'LineWidth',2); hold on;
plot(T,Ub_s,'--b','LineWidth',2); hold on;
plot(T,Ub_m,'--r','LineWidth',2); hold on;
plot(T,Ub_l,'--k','LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL','JS','JM','JL','QS','QM','QL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (m/s)')
%ylim([0 100])

%% speed based on mass & BL-Q10 in BL/s
figure(4)
plot(T,Ut_s/(1e-3*L_s),'b','LineWidth',2); hold on;
plot(T,Ut_m/(1e-3*L_m),'r','LineWidth',2); hold on;
plot(T,Ut_l/(1e-3*L_l),'k','LineWidth',2); hold on;
plot(T,Uw_s*ones(length(T),1)/(1e-3*L_s),'c','LineWidth',2); hold on;
plot(T,Uw_m*ones(length(T),1)/(1e-3*L_m),'m','LineWidth',2); hold on;
plot(T,Uw_l*ones(length(T),1)/(1e-3*L_l),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(T,Uj_s*ones(length(T),1)/(1e-3*L_s),'color',[0 0.5 0.75],'LineWidth',2); hold on;
plot(T,Uj_m*ones(length(T),1)/(1e-3*L_m),'color',[0.75 0 0.5],'LineWidth',2); hold on;
plot(T,Uj_l*ones(length(T),1)/(1e-3*L_l),'color',[0.33 0.33 0.33],'LineWidth',2); hold on;
plot(T,Ub_s/(1e-3*L_s),'--b','LineWidth',2); hold on;
plot(T,Ub_m/(1e-3*L_m),'--r','LineWidth',2); hold on;
plot(T,Ub_l/(1e-3*L_l),'--k','LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL','JS','JM','JL','QS','QM','QL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (BL/s)')

% Natural mortality per yr (Andersen & Beyer 2013)
%expoential coefficient
%Includes predation, excludes fishing
%mu = 0.35 * 4.5 * PI_s.^(-0.25);

