% Swimming speed
% Max recorded speeds are 56-110 km/h

L_s = 10^((log10(2)+log10(20))/2); 
L_m = 10^((log10(20)+log10(200))/2); 
L_l = 10^((log10(200)+log10(2000))/2); 

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

% Swimming speed 
%James' current (m/d)
%cm/s
U_s = (3.9*M_s.^0.13); %/100*60*60*24;
U_m = (3.9*M_m.^0.13); %/100*60*60*24;
U_l = (3.9*M_l.^0.13); %/100*60*60*24;

%Megrey et al with temp-dep (m/d)
T = -2:0.2:37;
T1= -2:0.2:8.9;
T2= 9:0.2:37;
Ut_s = (1.0*M_s.^0.0 * exp(0.0*T)); %/100*60*60*24;
Ut_m1 = (3.9*M_m.^0.13 * exp(0.149*T1)); %/100*60*60*24;
Ut_m2 = (15.0*M_m.^0.13 * exp(0.0*T2)); %/100*60*60*24;
Ut_l1 = (3.9*M_l.^0.13 * exp(0.149*T1)); %/100*60*60*24;
Ut_l2 = (15.0*M_l.^0.13 * exp(0.0*T2)); %/100*60*60*24;
Ut_m=[Ut_m1 Ut_m2];
Ut_l=[Ut_l1 Ut_l2];

%Q10 & BL (cm/s)
%Tref = 34.5490; %gives 0.1 BL/s at -2C, but never above 1.3BL/s
%Tref = 22.5617; %gives 3 BL/s at 40C, but <1 BL/s below 22.5
Tref=18.5; %median temp over historical period
Ub_s = L_s * exp(0.063*(T-Tref)) *0.1;
Ub_m = L_m * exp(0.063*(T-Tref)) *0.1;
Ub_l = L_l * exp(0.063*(T-Tref)) *0.1;

%%
figure(1)
plot(T,Ut_s,'b','LineWidth',2); hold on;
plot(T,Ut_m,'r','LineWidth',2); hold on;
plot(T,Ut_l,'k','LineWidth',2); hold on;
plot(T,U_s*ones(length(T),1),'c','LineWidth',2); hold on;
plot(T,U_m*ones(length(T),1),'m','LineWidth',2); hold on;
plot(T,U_l*ones(length(T),1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (cm/s)')

%%
figure(2)
plot(T,Ut_s/(0.1*L_s),'b','LineWidth',2); hold on;
plot(T,Ut_m/(0.1*L_m),'r','LineWidth',2); hold on;
plot(T,Ut_l/(0.1*L_l),'k','LineWidth',2); hold on;
plot(T,U_s*ones(length(T),1)/(0.1*L_s),'c','LineWidth',2); hold on;
plot(T,U_m*ones(length(T),1)/(0.1*L_m),'m','LineWidth',2); hold on;
plot(T,U_l*ones(length(T),1)/(0.1*L_l),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (BL/s)')

%%
figure(3)
plot(T,Ut_s,'b','LineWidth',2); hold on;
plot(T,Ut_m,'r','LineWidth',2); hold on;
plot(T,Ut_l,'k','LineWidth',2); hold on;
plot(T,U_s*ones(length(T),1),'c','LineWidth',2); hold on;
plot(T,U_m*ones(length(T),1),'m','LineWidth',2); hold on;
plot(T,U_l*ones(length(T),1),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(T,Ub_s,'--b','LineWidth',2); hold on;
plot(T,Ub_m,'--r','LineWidth',2); hold on;
plot(T,Ub_l,'--k','LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL','QS','QM','QL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (cm/s)')
ylim([0 100])

figure(4)
plot(T,Ut_s/(0.1*L_s),'b','LineWidth',2); hold on;
plot(T,Ut_m/(0.1*L_m),'r','LineWidth',2); hold on;
plot(T,Ut_l/(0.1*L_l),'k','LineWidth',2); hold on;
plot(T,U_s*ones(length(T),1)/(0.1*L_s),'c','LineWidth',2); hold on;
plot(T,U_m*ones(length(T),1)/(0.1*L_m),'m','LineWidth',2); hold on;
plot(T,U_l*ones(length(T),1)/(0.1*L_l),'color',[0.5 0.5 0.5],'LineWidth',2); hold on;
plot(T,Ub_s/(0.1*L_s),'--b','LineWidth',2); hold on;
plot(T,Ub_m/(0.1*L_m),'--r','LineWidth',2); hold on;
plot(T,Ub_l/(0.1*L_l),'--k','LineWidth',2); hold on;
legend('MS','MM','ML','WS','WM','WL','QS','QM','QL')
legend('location','northwest')
xlabel('Temp')
ylabel('Swimming speed (BL/s)')

% Natural mortality per yr (Andersen & Beyer 2013)
%expoential coefficient
%Includes predation, excludes fishing
%mu = 0.35 * 4.5 * PI_s.^(-0.25);

