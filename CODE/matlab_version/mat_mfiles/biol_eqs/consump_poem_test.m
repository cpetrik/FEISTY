%Test new consumption rates

clear all
close all

Lambda=0.7;

Sml_f.met = 0.03983707040056192;
Med_f.met = 0.017549368526009995;
Sml_p.met = 0.03983707040056192;
Med_p.met = 0.017549368526009995;
Lrg_p.met = 0.008662794702945285;
Sml_d.met = 0.03983707040056192;
Med_d.met = 0.017549368526009995;
Lrg_p.td = 0.6666666666666666;
Med_d.td = 0.5;
Sml_f.enc_zm = 0.056246807943889116;
Sml_p.enc_zm = 0.056246807943889116;
Sml_d.enc_d  = 1.1615590990689133e-6;
Med_f.enc_zl = 0.0013356065222719059;
Med_f.enc_zm = 0.0006903477046934407;
Med_p.enc_f  = 1.4256447382897658e-8;
Med_p.enc_p  = 1.4256447382897658e-8;
Med_d.enc_d  = 1.4256447382897658e-8;
Lrg_p.enc_f  = 4.6660571329057905e-11;
Lrg_p.enc_p  = 4.6660571329057905e-11;
Lrg_p.enc_d  = 2.3330285664528956e-11;
ENVR.Zm = 0.4842354382913072;
ENVR.Zl = 0.9368438618685104;
BENT.bio = 1.0e-5;
Sml_f.bio = 1.0e-5;
Sml_p.bio = 1.0e-5;
Sml_d.bio = 1.0e-5;
Med_d.bio = 1.0e-5;
Med_p.bio = 1.0e-5;
Med_f.bio = 1.0e-5;
Lrg_p.bio = 1.0e-5;
M_s=0.002529822128134705;
M_m=2.5298221281347053;
M_l=2529.822128134704;

%% consumption

Sml_f.con = Sml_f.enc_zm;
Sml_p.con = Sml_p.enc_zm;
Sml_d.con = Sml_d.enc_d; %+Sml_d.enc_zm;
Med_f.con = Med_f.enc_zm+Med_f.enc_zl;
Med_p.con = Med_p.enc_f+Med_p.enc_p;
Med_d.con = Med_d.enc_d;
Lrg_p.con = Lrg_p.enc_f+Lrg_p.enc_p+Lrg_p.enc_d;

%convert to rate by dividing predator biomass
Sml_f.con = Sml_f.con/Sml_f.bio;
Sml_p.con = Sml_p.con/Sml_p.bio;
Sml_d.con = Sml_d.con/Sml_d.bio;
Med_f.con = Med_f.con/Med_f.bio;
Med_p.con = Med_p.con/Med_f.bio;
Med_d.con = Med_d.con/Med_f.bio;
Lrg_p.con = Lrg_p.con/Lrg_p.bio;

%% nu = ((I/B)*Lambda) - met
Sml_f.nu = (Sml_f.con/Sml_f.bio)*Lambda - Sml_f.met;
Sml_p.nu = (Sml_p.con/Sml_p.bio)*Lambda - Sml_p.met;
Sml_d.nu = (Sml_d.con/Sml_d.bio)*Lambda - Sml_d.met;
Med_f.nu = (Med_f.con/Med_f.bio)*Lambda - Med_f.met;
Med_p.nu = (Med_p.con/Med_p.bio)*Lambda - Med_p.met;
Med_d.nu = (Med_d.con/Med_d.bio)*Lambda - Med_d.met;
Lrg_p.nu = (Lrg_p.con/Lrg_p.bio)*Lambda - Lrg_p.met;


%% nu = A*w^0.75;
%A = 10 gWW^0.25 per year at 15C)
A_s = 10%*M_s^0.25;
A_m = 10%*M_m^0.25;
A_l = 10%*M_l^0.25;

Sml_f.nu2 = A_s*M_s^0.75;
Sml_p.nu2 = A_s*M_s^0.75;
Sml_d.nu2 = A_s*M_s^0.75;
Med_f.nu2 = A_m*M_m^0.75;
Med_p.nu2 = A_m*M_m^0.75;
Med_d.nu2 = A_m*M_m^0.75;
Lrg_p.nu2 = A_l*M_l^0.75;

%%
n=3/4;
%h = 40; %g^(1-n)/yr at 10 degrees C
h = 40/365; %g^(1-n)/d at 10 degrees C
%Cmax = h*w^n;
Sml_f.Cmax = h*M_s^n;
Sml_p.Cmax = h*M_s^n;
Sml_d.Cmax = h*M_s^n;
Med_f.Cmax = h*M_m^n;
Med_p.Cmax = h*M_m^n;
Med_d.Cmax = h*M_m^n;
Lrg_p.Cmax = h*M_l^n;

%q = 0.8-1.0;
q=0.8;
%q=1.0;

%% Adjust gamma such that the feeding level f is around 0.6
gamma = 8000;
%V = gamma * w^q;
Sml_f.V = gamma * M_s^q;
Sml_p.V = gamma * M_s^q;
Sml_d.V = gamma * M_s^q;
Med_f.V = gamma * M_m^q;
Med_p.V = gamma * M_m^q;
Med_d.V = gamma * M_m^q;
Lrg_p.V = gamma * M_l^q;

%f = V*B / (Cmax + V*B);
Sml_f.f = Sml_f.V*ENVR.Zm / (Sml_f.Cmax + Sml_f.V*ENVR.Zm);
Sml_p.f = Sml_p.V*ENVR.Zm / (Sml_p.Cmax + Sml_p.V*ENVR.Zm);
Sml_d.f = Sml_d.V*ENVR.Zm / (Sml_d.Cmax + Sml_d.V*ENVR.Zm);
Med_f.f = Med_f.V*(ENVR.Zm+ENVR.Zl) / (Med_f.Cmax + Med_f.V*(ENVR.Zm+ENVR.Zl));
Med_p.f = Med_p.V*(Sml_f.bio+Sml_p.bio+Sml_d.bio) / (Med_p.Cmax + Med_p.V*(Sml_f.bio+Sml_p.bio+Sml_d.bio));
Med_d.f = Med_d.V*(Sml_f.bio+Sml_p.bio+Sml_d.bio) / (Med_d.Cmax + Med_d.V*(Sml_f.bio+Sml_p.bio+Sml_d.bio));
Lrg_p.f = Lrg_p.V*(Med_f.bio+Med_p.bio+Med_d.bio) / (Lrg_p.Cmax + Lrg_p.V*(Med_f.bio+Med_p.bio+Med_d.bio));

%Con = Cmax * f;
Sml_f.Con = Sml_f.Cmax * Sml_f.f;
Sml_p.Con = Sml_p.Cmax * Sml_p.f;
Sml_d.Con = Sml_d.Cmax * Sml_d.f;
Med_f.Con = Med_f.Cmax * Med_f.f;
Med_p.Con = Med_p.Cmax * Med_p.f;
Med_d.Con = Med_d.Cmax * Med_d.f;
Lrg_p.Con = Lrg_p.Cmax * Lrg_p.f;

%% handling time
%tau_s = 1/ (4*Sml_p.met);
tau_m = 1/ (4*Med_p.met);
tau_l = 1/ (4*Lrg_p.met);
tau_s = 1/ Sml_p.bio;

enc=0:1e-5:5e-4;
con = enc ./ (1+(tau_s.*enc));
con2 = Sml_p.bio * enc ./ (Sml_p.bio + enc);
figure
%subplot(2,2,1)
plot(enc,con,'b'); hold on;
plot(enc,con2,'--r')









