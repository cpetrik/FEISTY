% Visualize output of POEM
% Historic 1988-2010 with 3km model
% Time series plot of F vs. P with different Fishing rates

clear all
close all

% Fish data
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A020_Sm025_nmort1_BE08_noCC_RE00100';
%harv = 'pristine';
%tharv = 'F=0';
harv = 'All_fish01';
tharv = 'All fish F=0.1';
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
ppath = [pp cfile '/CC/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'Means_Historic3km_' harv '_' cfile '.mat']);


%% Plots in time
t = 1:length(sp_tmean); %time;
y = 1988 + (t-1)/12;

F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

%% Colors
cm4=[0 0 1;...    %b
    1 0 0;...     %r
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black

set(groot,'defaultAxesColorOrder',cm4);

%%
figure(1)
yyaxis left
plot(y,(P),'Linewidth',2);
ylabel('P Biomass (g m^-^2)')
xlim([y(1) y(end)])
yyaxis right
plot(y,(F),'Linewidth',2);
ylabel('F Biomass (g m^-^2)')
xlim([y(1) y(end)])
xlabel('Time (y)')
title(['Historic 3km ' tharv])
stamp(harv)
print('-dpng',[ppath 'Historic3km_',harv,'_FvP.png'])
 
