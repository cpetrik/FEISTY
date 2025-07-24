% Input data and params needed in advection-diffusion scheme

clear
close all

%% Velocities or Transports?
vpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CORE-forced/';
load([vpath,'Vel200_feb152013_run25_ocean_1988.mat'],'u','v');

% Grid
%1-D
load([vpath 'Data_grid_ocean_cobalt_ESM2Mcore.mat'],'GRD');
GRD1 = GRD;
clear GRD

%2-D
load([vpath 'Data_hindcast_grid_cp2D.mat'],'GRD')
GRD2 = GRD;
clear GRD

%% number of water cells
ID = GRD1.ID;
NX = length(ID);

% grid size
[ni,nj] = size(GRD2.mask);
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

mask = GRD2.mask;
dxtn = GRD2.dxtn;
dyte = GRD2.dyte;

%% load CORE spinup last mo for fish biomass
load('/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/CORE/Last_mo_Spinup1988_no_move_All_fish03_Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100.mat');

SFbio = nan*ones(ni,nj);
SPbio = nan*ones(ni,nj);
SDbio = nan*ones(ni,nj);
MFbio = nan*ones(ni,nj);
MPbio = nan*ones(ni,nj);
LPbio = nan*ones(ni,nj);

SFbio(ID) = Sml_f.bio;
SPbio(ID) = Sml_p.bio;
SDbio(ID) = Sml_d.bio;
MFbio(ID) = Med_f.bio;
MPbio(ID) = Med_p.bio;
LPbio(ID) = Lrg_p.bio;

%% load prey
load([vpath,'Data_ocean_cobalt_daily_1988.mat'],'COBALT');

%% define time
YEARS = 1;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
Mos = repmat(MNTH,1,YEARS);
tstep = 365 * 24 * 60 * 60; %time step in seconds

% Files to save
cname='CORE_1988_pelfish_dt1d_velDY';

MFbiov = zeros(NX,DAYS*YEARS);
LPbiov = zeros(NX,DAYS*YEARS);

%% call advec-diff
M=0;
n=0;
for Y=1:YEARS
    for mo = 1:length(Mos)
        M = M+1;
        for DAY = 1:Mos(mo)
            n=n+1;
            [num2str(mo) ',' num2str(DAY)];

            current(:,:,1) = u(:,:,DAY);
            current(:,:,2) = v(:,:,DAY);

            MZ = nan*ones(ni,nj);
            LZ = nan*ones(ni,nj);
            MZ(ID) = COBALT.Zm(:,DAY);
            LZ(ID) = COBALT.Zl(:,DAY);

            MFprey = SFbio + SFbio + SFbio + LZ + 0.1*MZ;
            LPprey = 0.5*MFbio + MPbio;

            MFbio = AdvectPredator( MFbio,...
                MFprey, ...
                current, ...
                tstep, ...
                dxtn, ...
                dyte, ...
                0.1, ...
                mask, ...
                nj, ...
                ni);

            LPbio = AdvectPredator( LPbio,...
                LPprey, ...
                current, ...
                tstep, ...
                dxtn, ...
                dyte, ...
                1.0, ...
                mask, ...
                nj, ...
                ni);

            MFbiov(:,n) = MFbio(ID);
            LPbiov(:,n) = LPbio(ID);
        end
    end
end

%% Save
spath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/CORE/';
save([spath 'AdvectPred_' cname '.mat'],...
    'MFbiov','LPbiov','GRD1','GRD2');

