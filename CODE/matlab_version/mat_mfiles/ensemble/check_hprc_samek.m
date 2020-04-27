% Test if laptop and HPRC results are the same

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

dp = ['/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];
hp = ['/Volumes/FEISTY/FEISTY_clim/Output/',...
    'Dc_cmax20-b250_D075_J100_A050_Sm025_nmort1_BE075_noCC_RE00100_Ka050/'];

%%
nfile = '/Volumes/FEISTY/NC/Matlab_new_size/param_ensemble/';
load([nfile 'LHS_param6_mid6_samek.mat'],'fx');

fsim = cell(length(fx),1);
sim = cell(length(fx),1);
for j = 1:length(fx)
    %! Change individual parameters
    pset = fx(j,:);
    %pset=round(pset,3);
    set_params6_samek(pset)
    
    %! Make core parameters/constants (global)
    const_params6_samek()
    
    %! Create a directory for output
    [fname,sname] = sub_fname_ensemble6_samek();
    fsim{j} = fname;
    sim{j} = sname;
end

%%
%same: 400, 476, 559, 640, 722
M = 722;
sfile = fsim{M};
sname = sim{M};
hfile = [hp fsim{M}(119:end)];

%% Last year means of laptop
load([sfile '.mat']);

[id,nt] = size(Clim_Bent.bio);
time=1:nt;
lyr=time((end-12+1):end);

CSPb=mean(Clim_Sml_p.bio(:,lyr),2);
CSFb=mean(Clim_Sml_f.bio(:,lyr),2);
CSDb=mean(Clim_Sml_d.bio(:,lyr),2);
CMPb=mean(Clim_Med_p.bio(:,lyr),2);
CMFb=mean(Clim_Med_f.bio(:,lyr),2);
CMDb=mean(Clim_Med_d.bio(:,lyr),2);
CLPb=mean(Clim_Lrg_p.bio(:,lyr),2);
CLDb=mean(Clim_Lrg_d.bio(:,lyr),2);
CBb=mean(Clim_Bent.bio(:,lyr),2);

CMFc=mean(Clim_Med_f.yield(:,lyr),2);
CMPc=mean(Clim_Med_p.yield(:,lyr),2);
CMDc=mean(Clim_Med_d.yield(:,lyr),2);
CLPc=mean(Clim_Lrg_p.yield(:,lyr),2);
CLDc=mean(Clim_Lrg_d.yield(:,lyr),2);

%% Last year means of HPRC
load([hfile '.mat']);
[id,nt] = size(Clim_Bent.bio);
time=1:nt;
lyr=time((end-12+1):end);

HSPb=mean(Clim_Sml_p.bio(:,lyr),2);
HSFb=mean(Clim_Sml_f.bio(:,lyr),2);
HSDb=mean(Clim_Sml_d.bio(:,lyr),2);
HMPb=mean(Clim_Med_p.bio(:,lyr),2);
HMFb=mean(Clim_Med_f.bio(:,lyr),2);
HMDb=mean(Clim_Med_d.bio(:,lyr),2);
HLPb=mean(Clim_Lrg_p.bio(:,lyr),2);
HLDb=mean(Clim_Lrg_d.bio(:,lyr),2);
HBb=mean(Clim_Bent.bio(:,lyr),2);

HMFc=mean(Clim_Med_f.yield(:,lyr),2);
HMPc=mean(Clim_Med_p.yield(:,lyr),2);
HMDc=mean(Clim_Med_d.yield(:,lyr),2);
HLPc=mean(Clim_Lrg_p.yield(:,lyr),2);
HLDc=mean(Clim_Lrg_d.yield(:,lyr),2);

%% Compare
sf = sum(CSFb==HSFb)
sp = sum(CSPb==HSPb)
sd = sum(CSDb==HSDb)
mf = sum(CMFb==HMFb)
mp = sum(CMPb==HMPb)
md = sum(CMDb==HMDb)
lp = sum(CLPb==HLPb)
ld = sum(CLDb==HLDb)
b  = sum(CBb ==HBb)

%%
mfc = sum(CMFc==HMFc)
mpc = sum(CMPc==HMPc)
mdc = sum(CMDc==HMDc)
lpc = sum(CLPc==HLPc)
ldc = sum(CLDc==HLDc)
