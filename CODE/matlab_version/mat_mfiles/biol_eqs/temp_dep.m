% Temperature-dependence

clear all
close all

%% Megrey et al.
T = -2:30;

% Consumption  -------------------------
% age-0
te1 = 1; 
te2 = 15;
te3 = 17;
te4 = 23;
xk1 = 0.1;
xk2 = 0.98;
xk3 = 0.98;
xk4 = 0.01;

tt5 = 1 ./ (te2 - te1);
t5 = tt5 .* log((xk2.*(1-xk1)) ./ (xk1.*(1-xk2)));
t4 = exp(t5.*(T-te1));
tt7 = 1 ./ (te4-te3);
t7 = tt7 .* log((xk3.*(1-xk4)) ./ (xk4.*(1-xk3)));
t6 = exp(t7.*(te4-T));
gcta0 = (xk1.*t4) ./ (1+xk1.*(t4-1));
gctb0 = (xk4.*t6) ./ (1+xk4.*(t6-1));
fT0 = gcta0.*gctb0;
%Cmax = a.*W^(-b).*fT;

% age-1
te1 = 1; 
te2 = 15;
te3 = 17;
te4 = 25;
xk1 = 0.1;
xk2 = 0.98;
xk3 = 0.98;
xk4 = 0.01;

tt5 = 1 ./ (te2 - te1);
t5 = tt5 .* log((xk2.*(1-xk1)) ./ (xk1.*(1-xk2)));
t4 = exp(t5.*(T-te1));
tt7 = 1 ./ (te4-te3);
t7 = tt7 .* log((xk3.*(1-xk4)) ./ (xk4.*(1-xk3)));
t6 = exp(t7.*(te4-T));
gcta1 = (xk1.*t4) ./ (1+xk1.*(t4-1));
gctb1 = (xk4.*t6) ./ (1+xk4.*(t6-1));
fT1 = gcta1.*gctb1;
%Cmax = a.*W^(-b).*fT;

% age-2
te1 = 1; 
te2 = 13;
te3 = 15;
te4 = 23;
xk1 = 0.1;
xk2 = 0.98;
xk3 = 0.98;
xk4 = 0.01;

tt5 = 1 ./ (te2 - te1);
t5 = tt5 .* log((xk2.*(1-xk1)) ./ (xk1.*(1-xk2)));
t4 = exp(t5.*(T-te1));
tt7 = 1 ./ (te4-te3);
t7 = tt7 .* log((xk3.*(1-xk4)) ./ (xk4.*(1-xk3)));
t6 = exp(t7.*(te4-T));
gcta2 = (xk1.*t4) ./ (1+xk1.*(t4-1));
gctb2 = (xk4.*t6) ./ (1+xk4.*(t6-1));
fT2 = gcta2.*gctb2;
%Cmax = a.*W^(-b).*fT;

% plot
figure(1)
plot(T,fT0,'r')
hold on
plot(T,fT1,'b')
hold on
plot(T,fT2,'k')

% Whole range
fT = gcta2.*gctb1;
figure(2)
plot(T,fT,'m')

% Metabolism/Respiration -------------------------
c0 = 0.083;
c12 = 0.0548;

fR0 = exp(c0*T);
fR12 = exp(c12*T);
%R = a.*W^(-b).*fR .* activity .*5.258;

figure(3)
plot(T,fR0,'r')
hold on
plot(T,fR12,'k')


%% Rall et al 2012 Consumption
%k = 1.38*1e-23 * (1/1.60e-19); %eV K^-1
k = 8.617e-5; %eV K^-1
T0 = 293.15;
a0 = -18.77;
bi = 0.93;
bj = 0.15;
Eh = 0.40;
h0 = 4.54;
ci = -0.02;
cj = -0.14;
Ea = 1.15;

% Varying temp
Tkel = T + 273.25;
mi = 0.0372;
mj = 2.53e3;
laT = a0 + bi.*log(mi) + bj.*log(mj) + Ea.*(Tkel-T0)./(k.*Tkel.*T0);
lhT = h0 + ci.*log(mi) + cj.*log(mj) + Eh.*(Tkel-T0)./(k.*Tkel.*T0);
aT = exp(laT);
hT = exp(lhT);

figure(4)
subplot(2,1,1)
plot(T,aT);
ylabel('attack rate')
subplot(2,1,2)
plot(T,hT);
ylabel('handling time')
xlabel('temperature')

% Varying prey size
Tkel = 15 + 273.25;
mi = linspace(1.5e-4,9.3,20);
mj = 2.53e3;
laP = a0 + bi.*log(mi) + bj.*log(mj) + Ea.*(Tkel-T0)./(k.*Tkel.*T0);
lhP = h0 + ci.*log(mi) + cj.*log(mj) + Eh.*(Tkel-T0)./(k.*Tkel.*T0);
aP = exp(laP);
hP = exp(lhP);

figure(5)
subplot(2,1,1)
plot(log(mi),aP);
ylabel('attack rate')
subplot(2,1,2)
plot(log(mi),hP);
ylabel('handling time')
xlabel('ln prey weight (mg)')

% Varying fish size
Tkel = 15 + 273.25;
mi = 0.0372;
mj = linspace(0.8,8e7,20);
laF = a0 + bi.*log(mi) + bj.*log(mj) + Ea.*(Tkel-T0)./(k.*Tkel.*T0);
lhF = h0 + ci.*log(mi) + cj.*log(mj) + Eh.*(Tkel-T0)./(k.*Tkel.*T0);
aF = exp(laF);
hF = exp(lhF);

figure(6)
subplot(2,1,1)
plot(log(mj),aF);
ylabel('attack rate')
subplot(2,1,2)
plot(log(mj),hF);
ylabel('handling time')
xlabel('ln fish weight (mg)')

% Varying prey abundance
Tkel = 15 + 273.25;
mi = 0.0372;
mj = 2.53e3;
%Mesozoop biom ranges from 1e2 to 1e4 mg C/m2
% 1e4 mg C = 9e4 mg wet wgt
% smallest zoop of 1.5e-4 mg -> abund = 6.1e8 per m2
N = linspace(0,1e9,100);
la = a0 + bi.*log(mi) + bj.*log(mj) + Ea.*(Tkel-T0)./(k.*Tkel.*T0);
lh = h0 + ci.*log(mi) + cj.*log(mj) + Eh.*(Tkel-T0)./(k.*Tkel.*T0);
a = exp(la);
h = exp(lh);
Fmax = 1/h;
f = a.*N ./ (Fmax + a.*N);
C = Fmax .* f;

figure(7)
plot(N,C)
ylabel('consumption (s^-^1)')
xlabel('prey abundance')

figure(8)
plot(N,C*60*60*24)
ylabel('consumption (d^-^1)')
xlabel('prey abundance')


%% Pope et al 2009 Temp correction factor
E = 0.6;
k = 8.617343e-5;
T = -2:30;
Tkel = T + 273.25;
tau1 = 5 + 273.25;
tau2 = 15 + 273.25;
tau3 = 25 + 273.25;
omg1 = exp(-(E./k).*(1./Tkel - 1./tau1));
omg2 = exp(-(E./k).*(1./Tkel - 1./tau2));
omg3 = exp(-(E./k).*(1./Tkel - 1./tau3));

figure(9)
plot(T,omg1,'b'); hold on;
plot(T,omg2,'k'); hold on;
plot(T,omg3,'r'); hold on;
legend('5C','15C','25C')
xlabel('Temp')
ylabel('correction factor')







