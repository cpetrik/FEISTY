
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

q=0.8;
h=85;
w=0.0025298;
prey=0:0.005:0.3;
gamma=8e4;
% Small
AS = (gamma*M_s^0.8)/365;
cmaxS = (h.*M_s^0.75)/365 / M_s;
metS = (0.05 * cmaxS);
%to get con/cmax = 0.6
lam1S = 365*1.5*(cmaxS) ./ (prey.*M_s^q);
figure
plot(prey,lam1S,'LineWidth',2)
%to get con*assim>metab
lam2S = 365* (metS) ./ ((0.8-0.05) .* M_s^q .*prey);
figure
plot(prey,lam2S)

%%
AS1=(1e6*M_s^0.8)/365;
AS2=(4e6*M_s^0.8)/365;
AS3=(8e6*M_s^0.8)/365;
con1=(cmaxS.*AS1*prey)./(cmaxS+AS1.*prey);
con2=(cmaxS.*AS2*prey)./(cmaxS+AS2.*prey);
con3=(cmaxS.*AS3*prey)./(cmaxS+AS3.*prey);
figure
plot(prey,con1,'b'); hold on;
plot(prey,con2,'k'); hold on;
plot(prey,con3,'r'); hold on;


%% Med
AM = (gamma*M_m^0.8)/365;
cmaxM = (h.*M_m^0.75)/365 / M_m;
metM = (0.05 * cmaxM);
%to get con/cmax = 0.6
lam1M = 365*1.5*(cmaxM) ./ (prey.*M_m^q);
figure
plot(prey,lam1M,'LineWidth',2)
%to get con*assim>metab
lam2M = 365* (metM) ./ ((0.8-0.05) .* M_m^q .*prey);
figure
plot(prey,lam2M)

% Lrg
AL = (gamma*M_l^0.8)/365;
cmaxL = (h.*M_l^0.75)/365 / M_l;
metL = (0.05 * cmaxL);
%to get con/cmax = 0.6
lam1L = 365*1.5*(cmaxL) ./ (prey.*M_l^q);
figure
plot(prey,lam1L,'LineWidth',2)
%to get con*assim>metab
lam2L = 365* (metL) ./ ((0.8-0.05) .* M_l^q .*prey);
figure
plot(prey,lam2L)

%% Temp effects
temp=-2:30;
cmaxs = exp(0.063*(temp-10.0)) .* (h.*M_s^0.75)/365 / M_s;
cmaxm = exp(0.063*(temp-10.0)) .* (h.*M_m^0.75)/365 / M_m;
cmaxl = exp(0.063*(temp-10.0)) .* (h.*M_l^0.75)/365 / M_l;

figure
subplot(3,1,1)
plot(temp,cmaxs)
subplot(3,1,2)
plot(temp,cmaxm)
subplot(3,1,3)
plot(temp,cmaxl)



