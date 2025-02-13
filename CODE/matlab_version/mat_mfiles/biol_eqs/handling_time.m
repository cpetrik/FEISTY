clear all
close all

%!Individual Mass (g)
M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);
m=[M_s; M_m; M_l];

temp=15;
Ccmax = (exp(0.063*(temp-10.0)) .* 20 .* m.^(-0.25)) ./365.0;

%% handling time
% Sml_f.met = 0.03983707040056192;
% Med_f.met = 0.017549368526009995;
% Sml_p.met = 0.03983707040056192;
% Med_p.met = 0.017549368526009995;
% Lrg_p.met = 0.008662794702945285;

% tau_s = 1/ (Sml_p.met);
% tau_m = 1/ (Med_p.met);
% tau_l = 1/ (Lrg_p.met);

tau = 1./Ccmax;

enc=0:0.1:7;
% con = enc ./ (1+(tau_s.*enc));
con = enc ./ (1+(tau.*enc));

figure(1)
%subplot(2,2,1)
plot(enc,con(1,:),'r'); hold on;
plot(enc,con(2,:),'k'); hold on;
plot(enc,con(3,:),'b'); hold on;

%% Mediums

con0 = Ccmax(2)*enc ./ (Ccmax(2)+enc);
con1 = 0.1*Ccmax(2)*enc ./ (0.1*Ccmax(2)+enc);
con2 = 10*Ccmax(2)*enc ./ (10*Ccmax(2)+enc);
con3 = 0.1*Ccmax(2)*enc ./ (Ccmax(2)+0.1*enc);
con4 = 10*Ccmax(2)*enc ./ (Ccmax(2)+10*enc);

figure(2)
%subplot(2,2,1)
plot(enc,con0,'k'); hold on;
plot(enc,con1,'b'); hold on;
plot(enc,con2,'r'); hold on;

figure(3)
%subplot(2,2,1)
plot(enc,con0,'k'); hold on;
plot(enc,con3,'b'); hold on;
plot(enc,con4,'r'); hold on;
xlim([0 1])
%%
con1 = 0.25*Ccmax(2)*enc ./ (Ccmax(2)+0.25*enc);
con2 = 0.50*Ccmax(2)*enc ./ (Ccmax(2)+0.50*enc);
con3 = 0.75*Ccmax(2)*enc ./ (Ccmax(2)+0.75*enc);
con4 = 0.90*Ccmax(2)*enc ./ (Ccmax(2)+0.90*enc);

figure(4)
%subplot(2,2,1)
plot(enc,con0,'k'); hold on;
plot(enc,con1,'b'); hold on;
plot(enc,con2,'r'); hold on;
plot(enc,con3,'c'); hold on;
plot(enc,con4,'m'); hold on;
%xlim([0 1])

