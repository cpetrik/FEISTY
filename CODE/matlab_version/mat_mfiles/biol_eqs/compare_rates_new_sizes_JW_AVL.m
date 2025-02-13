% Compare diff Cmax, metab, A, ingest eqs

clear all
close all

%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

%! Body lengths (mm)
% Convert from mm to cm and use their const coeff = 0.01g/cm3
L_s = 10.0 * (M_s/0.01)^(1/3); % small
L_m = 10.0 * (M_m/0.01)^(1/3); % medium
L_l = 10.0 * (M_l/0.01)^(1/3); % large

m=[M_s; M_m; M_l];
len=[L_s; L_m; L_l];

%Hartvig et al. constants
q=1-0.8;
h=85;
w=0.0025298;
gamma=0.8e4;
n=1-0.75;
p=1-0.75;
k=10;
u=0.84;

%mizer constants
h2=40;
gamma2=2.9e3;
k2=4.8;
u2=3;

%J&C 15 constants
Q=1-0.9;
H=25;
W=0.0025298;
G=1.37e4;
N=1-0.67;
K=2.5;
U=0.5;

prey=0:0.005:0.3;
temp=15;

%% Mortality
%Hartvig
Hu = exp(0.063*(temp-10.0)) .* (u.*m.^(-n))/365;

%mizer
Zu = exp(0.063*(temp-10.0)) .* (u2.*m.^(-n))/365;

%J&C
Ju = exp(0.063*(temp-10.0)) .* (U.*m.^(-N))/365;

%JW
Wu = [0; 0.05; 0.1] /365;
    
%AVL

%% Swim speed
U = ((3.9*m.^0.13)/100*60*60*24);
Ut = ((3.9*m.^0.13 * exp(0.149*temp)) /100*60*60*24);
	
%% Respiration
%Hartvig
Hmet = exp(0.063*(temp-10.0)) .* k.*m.^-p;

%mizer
Zmet = exp(0.063*(temp-10.0)) .* k2.*m.^-p;

%J&C
Jmet = exp(0.063*(temp-10.0)) .* K.*m.^-p;

%Me
Mcmax = exp(0.063*(temp-15.0)) .* (60.*m.^-n)/365;
Mmet = exp(0.063*(temp-15.0)) .* 0.4*Mcmax;

%JW
%Resting
WmetP = 0.0033*m.^-0.13;
WmetD = 0.00033*m.^-0.13;
%Active
WmetA = exp(0.03*(3.9*m.^0.13));
WmetAP = WmetP * 5.258;
WmetAD = WmetD * 5.258;
WmetAP2 = WmetP .* exp(0.03*(U*100/60/60/24)) * 5.258;
WmetAD2 = WmetD .* exp(0.03*(U*100/60/60/24)) * 5.258;
WmetAP2t = WmetP .* exp(0.03*(Ut*100/60/60/24)) * 5.258;
WmetAD2t = WmetD .* exp(0.03*(Ut*100/60/60/24)) * 5.258;
    
%AVL

%% Max ingestion
%Kiorboe & Hirst
Imax = exp(0.063*(temp-15.0)) .* 2.5 .* m.^(-0.51) .* 24e-3;

%Hartvig
Hcmax = exp(0.063*(temp-10.0)) .* (h.*m.^-n)/365;

%mizer
Zcmax = exp(0.063*(temp-10.0)) .* (h2.*m.^-n)/365;

%J&C
Jcmax = exp(0.063*(temp-10.0)) .* (H.*m.^-N)/365;

%Me
Mcmax = exp(0.063*(temp-15.0)) .* (60.*m.^-n)/365;

%JW
WcmaxP = (4*WmetP);
WcmaxD = (4*WmetD);
WcmaxA = (4*WmetA);
WcmaxAP = (4*WmetAP);
WcmaxAD = (4*WmetAD);
WcmaxAP2 = (4*WmetAP2);
WcmaxAD2 = (4*WmetAD2);
WcmaxAP2t = (4*WmetAP2t);
WcmaxAD2t = (4*WmetAD2t);

%AVL


%% Max clearance rate/search volume
%Kiorboe & Hirst = me
Fmax = exp(0.063*(temp-15.0)) .* 10^(3.24) .* m.^(-0.24) .* (24e-3/9);

%Hartvig
HA = exp(0.063*(temp-10.0)) .* (gamma*m.^-q)/365;

%mizer
ZA = exp(0.063*(temp-10.0)) .* (gamma2*m.^-q)/365;

%J&C
JA = exp(0.063*(temp-10.0)) .* (G*m.^-Q)/365;

%JW
Tu = [1.0;0.5;0.1];
WA = (U .* ((len./1000)*3) .* Tu) ./ m;
WAt = (Ut .* ((len./1000)*3) .* Tu) ./ m;

%AVL

%%
figure(1)
%subplot(2,2,1)
plot(log10(m),log10(Fmax),'r','MarkerSize',25); hold on;
plot(log10(m),log10(HA),'b','MarkerSize',25); hold on;
plot(log10(m),log10(JA),'k','MarkerSize',25); hold on;
plot(log10(m),log10(ZA),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(WA),'c','MarkerSize',25); hold on;
plot(log10(m),log10(WAt),'g','MarkerSize',25); hold on;
ylabel('log10 clearance (m^3 gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','Hart','J&C','mizer','JW','JWt')
legend('location','southwest')
title('Specific clearance rate')
print('-dpng','Clearance_new_sizes_JW.png')

figure(2)
%subplot(2,2,1)
plot(log10(m),log10(Imax),'r','MarkerSize',25); hold on;
plot(log10(m),log10(Hcmax),'b','MarkerSize',25); hold on;
plot(log10(m),log10(Jcmax),'k','MarkerSize',25); hold on;
plot(log10(m),log10(Zcmax),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mcmax),'ob','MarkerSize',10); hold on;
plot(log10(m),log10(WcmaxP),'color',[0 1 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxD),'color',[0 1 0],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxA),'color',[1 1 0],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxAP),'color',[1 0.5 0],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxAD),'color',[1 0 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxAP2),'color',[0.5 0 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxAD2),'color',[0 0.5 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxAP2t),'color',[0.25 0.25 0.25],'MarkerSize',25); hold on;
plot(log10(m),log10(WcmaxAD2t),'color',[0.75 0.75 0.75],'MarkerSize',25); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('K&H','Hart','J&C','mizer','me','P','D','A','AP','AD','AP2','AD2',...
    'AP2t','AD2t')
legend('location','eastoutside')
title('Specific ingestion rate')
print('-dpng','Ingestion_new_sizes_JW.png')

figure(3)
%subplot(2,2,1)
plot(log10(m),log10(Hmet),'b','MarkerSize',25); hold on;
plot(log10(m),log10(Jmet),'k','MarkerSize',25); hold on;
plot(log10(m),log10(Zmet),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),log10(Mmet),'ob','MarkerSize',10); hold on;
plot(log10(m),log10(WmetP),'color',[0 1 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetD),'color',[0 1 0],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetA),'color',[1 1 0],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetAP),'color',[1 0.5 0],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetAD),'color',[1 0 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetAP2),'color',[0.5 0 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetAD2),'color',[0 0.5 1],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetAP2t),'color',[0.25 0.25 0.25],'MarkerSize',25); hold on;
plot(log10(m),log10(WmetAD2t),'color',[0.75 0.75 0.75],'MarkerSize',25); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','me','P','D','A','AP','AD','AP2','AD2',...
    'AP2t','AD2t')
legend('location','eastoutside')
title('Specific respiration rate')
print('-dpng','Respiration_new_sizes_JW.png')

figure(4)
plot(log10(m),(Hu),'b','MarkerSize',25); hold on;
plot(log10(m),(Ju),'k','MarkerSize',25); hold on;
plot(log10(m),(Zu),'color',[0.5 0.5 0.5],'MarkerSize',25); hold on;
plot(log10(m),(Wu),'c','MarkerSize',25); hold on;
ylabel('mortality rate (d^-^1)')
xlabel('log10 weight (gWW)')
legend('Hart','J&C','mizer','JW')
legend('location','southwest')
title('Daily mortality rate')
print('-dpng','Mortality_new_sizes_JW.png')


