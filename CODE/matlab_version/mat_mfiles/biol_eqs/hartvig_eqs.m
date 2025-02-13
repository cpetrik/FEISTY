% Hartvig et al eqs

clear all
close all

%% Size & temp
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

%temp
T = -2:0.2:37;
%phi = prey biomass available
phi = 0.002:0.01:2;

%% Constants
q = 0.8;
y = 0.8 * 10^4; %m^3 g^(?q)/year at 10C
%y = 4.3 * 10^4; %m^3 g^(?q)/year at 10C
h = 85; %g^(1-n)/yr at 10C
n = 0.75;
a = 0.6;
k = 10; %g^(1-p)/yr at 10C
p = 0.75;

%% Biol
m = M_s;

%Search volume:
v = y*m^q; %m^3

%Encountered food:
enc = v*phi;
figure(1)
plot(phi,enc)
title('Encounter per individual')
ylabel('Encountered biomass per individual (g/yr)')
xlabel('Prey concentration (g/m^2)')

%Feeding level:
f1 = (v.*phi) ./ (v.*phi + h.*m^n);
f2 = enc ./ (enc+h.*m^n);
figure(2)
plot(phi,f1,'k'); hold on;
plot(phi,f2,'b'); hold on;
ylabel('Feeding level')
xlabel('Prey concentration (g/m^2)')

%Energy available for growth:
E1 = a.*f1.*h.*m^n - k.*m^p;
E2 = a.*f2.*h.*m^n - k.*m^p;
%E3 = ((a.*cmax.*enc) ./ (enc+cmax)) .* m^n - k.*m^p;
figure(3)
plot(phi,E1,'k'); hold on;
plot(phi,E2,'b'); hold on;
ylabel('Energy available for growth (g)')
xlabel('Prey concentration (g/m^2)')

%Standard and active metabolism:
metab = k*m^p;

%Critical feeding level:
fc1 = (k.*m^(n-p)) ./ (a.*h);
fc2 = (k./(a.*h));
%fc ~ 0.2

%% POEM code vs. prey
tprey=1;
fcrit=0.3;
temp=10;
wgt=M_s;
prey=phi;
pred=M_s;

%Cmax
cmax = (exp(0.063*(temp-10.0)) * h * wgt^(3/4)) /365; %h value for temp=10C

%Metabolism
bas = (fcrit * cmax);
act = (exp(0.063*(temp-10.0)) * k * wgt^(3/4)) /365; %Charlie thinks another way is relate it to amount consumed
met = bas + act;
bvec = bas*ones(size(prey));
avec = met*ones(size(prey));

%Encounter rates
A = (exp(0.063*(temp-10.0)) * y * wgt^(q)) /365;   %coeffs for per yr -> per day
enc = pred*prey*A*tprey;
con = cmax .* enc ./ (cmax + enc);

figure(4)
subplot(2,2,1)
plot(prey,enc)
title('Encounter per larva')
ylabel('Encountered biomass (g/d)')
xlabel('Zoop concentration (g/m^2)')

figure(5)
subplot(2,2,1)
plot(prey,con,'k'); hold on;
plot(prey,bvec,'b'); hold on;
plot(prey,avec,'r'); hold on;
legend('con','bas','act')
title('Consumption per larva')
ylabel('Consumed biomass (g/d)')
xlabel('Zoop concentration (g/m^2)')

figure(6)
subplot(2,2,1)
plot(enc,con,'k'); hold on;
plot(enc,bvec,'b'); hold on;
plot(enc,avec,'r'); hold on;
legend('con','bas','act')
title('Consumption per larva')
ylabel('Consumed biomass (g/d)')
xlabel('Encountered biomass (g/m^2)')

%% POEM code vs Temp
tprey=1;
fcrit=0.3;
temp=T;
wgt=M_s;
prey1=0.02;
prey2=0.2;
prey3=2.0;
pred=M_s;
y=0.8e4;

%Cmax
cmax = (exp(0.063.*(temp-10.0)) .* h .* wgt^(3/4)) ./365; %h value for temp=10C

%Metabolism
bas = (fcrit .* cmax);
act = (exp(0.063.*(temp-10.0)) .* k .* wgt^(3/4)) ./365; %Charlie thinks another way is relate it to amount consumed
met = bas + act;

%Encounter rates
A = (exp(0.063.*(temp-10.0)) .* y .* wgt^(q)) ./365;   %coeffs for per yr -> per day
enc1 = pred.*prey1.*A.*tprey;
enc2 = pred.*prey2.*A.*tprey;
enc3 = pred.*prey3.*A.*tprey;
con1 = cmax .* enc1 ./ (cmax + enc1);
con2 = cmax .* enc2 ./ (cmax + enc2);
con3 = cmax .* enc3 ./ (cmax + enc3);

figure(7)
subplot(2,2,1)
plot(temp,enc1,'b'); hold on;
plot(temp,enc2,'k'); hold on;
plot(temp,enc3,'r'); hold on;
legend('min','mean','max')
title('Encounter per larva')
ylabel('Encountered biomass (g/d)')
xlabel('Temp (C)')

figure(8)
subplot(2,2,1)
plot(temp,con1,'b'); hold on;
plot(temp,con2,'k'); hold on;
plot(temp,con3,'r'); hold on;
plot(temp,bas,'c'); hold on;
plot(temp,met,'m'); hold on;
legend('min','mean','max','bas','act')
title('Consumption per larva')
ylabel('Consumed biomass (g/d)')
xlabel('Temp (C)')

%%
f0=0.6;
beta=100;
lam=2-n+q;
theta=0:0.1:1;

gam = (f0 * h * beta^(2-lam)) ./ ((1-f0) * theta * sqrt(2*pi));

figure(9)
subplot(2,2,1)
plot(theta,gam,'b'); hold on;
