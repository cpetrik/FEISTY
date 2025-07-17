% Test 1st order upwind, 2nd order centered, 4th order centered 
% advection schemes from MOM5
% Input data and params needed in advection-diffusion scheme

clear all
close all

% Transports path
vpath = '/Volumes/GFDL/POEM_JLD/esm2m_hist/';

% Grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_hindcast_grid_cp2D.mat')

%% number of water cells
ID = find(GRD.mask==1);
NX = length(ID);

% grid size
[ni,nj] = size(GRD.mask);
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

% define time
YEARS = 1;
DAYS = 365;
tstep = 1; %time step in hours

tdt = num2str(1000+int64(100*tstep));

%% define a patch to advect
bio = zeros(ni,nj);
%File name to save

%Global
% bio = 100*ones(ni,nj);   %Global
% cname='Global_even_dt1hr_esm2m2000_vel_b100_area';

%Atl-Arctic
bio(220:240,:) = 1.0e2; bio(121:141,195:200) = 1.0e2; 
cname=['AtlArc_even_dt',tdt(2:end),'hr_esm2m2000_diffv2'];

%seed equator
% bio(:,84:109) = 1.0e2;     
% cname='Eq_even_dt1hr_esm2m2000_vel_b100_area';

%seed Atl
%bio(220:240,:) = 1.0e2;    
%cname='Atl_even_dt1hr_esm2m2000_vel_b100_area';

%seed Pac
%bio(59:79,:) = 1.0e2;      
%cname='Pac_even_dt1hr_esm2m2000_vel_b100_area';

%seed Indian W
%bio(5:25,:) = 1.0e2;      
%cname='WInd_even_dt1hr_esm2m2000_vel_b100_area';

%seed Indian E
%bio(340:360,:) = 1.0e2;    
%cname='EInd_even_dt1hr_esm2m2000_vel_b100_area';

%seed Arctic
%bio(:,181:200) = 1.0e2;    
%cname='Arc_even_dt1hr_esm2m2000_vel_b100_area';

%seed Antarctic
%bio(:,12:32) = 1.0e2;      
%cname='Ant_even_dt1hr_esm2m2000_vel_b100_area';

%w/i Natasha Atl
%bio(206:295,150:177) = 1.0e2; 
%cname='NEAtl_even_dt1hr_esm2m2000_vel_b100_area';

bio = bio .* GRD.mask;

biov = zeros(NX,DAYS*YEARS);

%% do advec-diff
% % define diffusivity
% %horizontal diffusivity varies by two orders of magnitude from 10^2 to 10^4 m-2 s-1 
% K = 600.0;
% 
% bio1 = bio;
% bio2 = bio;
% bio4 = bio;
% biov1 = zeros(NX,DAYS*YEARS);
% biov2 = zeros(NX,DAYS*YEARS);
% biov4 = zeros(NX,DAYS*YEARS);
% 
% n=0;
% for Y=1:YEARS
%     yr = num2str(Y+1988-1);
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         DAY
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         %1st order upwind advect
%         bio1 = sub_advec1_diff_vel(GRD,bio1,K,U,V,ni,nj,tstep);
%         biov1(:,n) = bio1(ID);
%         %2nd order centered advect
% %         bio2 = sub_advec2_diff_vel(GRD,bio2,K,U,V,ni,nj,tstep);
% %         biov2(:,n) = bio2(ID);
% %         %4th order centered advect
% %         bio4 = sub_advec4_diff_vel(GRD,bio4,K,U,V,ni,nj,tstep);
% %         biov4(:,n) = bio4(ID);
%     end
% end
% 
% % Save
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_diff_1st_2nd_4th_K',num2str(K),...
% '_' cname '.mat'],'biov1','biov2','biov4','biov','bio2','bio4');

%% do advec only
% % 1ST ORDER
% bio1 = bio;
% biov1 = zeros(NX,DAYS*YEARS);
% n=0;
% for Y=1:YEARS
%     %yr = num2str(Y+1988-1);
%     yr = num2str(Y+1996-1); %only have 1996-2000
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         [num2str(DAY) ',1']
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         bio1 = sub_advec_vel_1st_order(GRD,bio1,U,V,ni,nj,tstep);
%         biov1(:,n) = bio1(ID);
%     end
% end
% 
% % 2ND ORDER
% bio2 = bio;
% biov2 = zeros(NX,DAYS*YEARS);
% n=0;
% for Y=1:YEARS
%     %yr = num2str(Y+1988-1);
%     yr = num2str(Y+1996-1); %only have 1996-2000
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         [num2str(DAY) ',2']
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         bio2 = sub_advect_vel_2nd_order(GRD,bio2,U,V,ni,nj,tstep);
%         biov2(:,n) = bio2(ID);
%     end
% end
% 
% % 4TH ORDER
% bio4 = bio;
% biov4 = zeros(NX,DAYS*YEARS);
% n=0;
% for Y=1:YEARS
%     %yr = num2str(Y+1988-1);
%     yr = num2str(Y+1996-1); %only have 1996-2000
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         [num2str(DAY) ',4']
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         bio4 = sub_advect_vel_4th_order(GRD,bio4,U,V,ni,nj,tstep);
%         biov4(:,n) = bio4(ID);
%     end
% end
% 
% % Save
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_1st_2nd_4th_' cname '.mat'],...
%     'biov1','biov2','biov4','biov','bio2','bio4');

%% do diff only
%Global
% bio = 100*rand(ni,nj);   
% cname='Global_rand_dt1hr_esm2m2000_vel_b100_area';
%Gradient highest at eq
% g1 = 1:100;
% g2 = 100:-1:1;
% bio = zeros(ni,nj);
% bio(:,1:100) = repmat(g1,ni,1);
% bio(:,101:200) = repmat(g2,ni,1);
% cname='Eq_grad_dt1hr_esm2m2000_vel_b100_area';
%Gradient highest at poles
% g1 = 1:100;
% g2 = 100:-1:1;
% bio = zeros(ni,nj);
% bio(:,1:100) = repmat(g2,ni,1);
% bio(:,101:200) = repmat(g1,ni,1);
% cname='Pole_grad_dt1hr_esm2m2000_vel_b100_area';

% bio = bio .* GRD.mask;

bio1 = bio;
bio2 = bio;
bio4 = bio;
biov1 = zeros(NX,DAYS*YEARS);
biov2 = zeros(NX,DAYS*YEARS);
biov4 = zeros(NX,DAYS*YEARS);

% define diffusivity
K = 600.0;
U = zeros(ni,nj);
V = zeros(ni,nj);
n=0;
for Y=1:YEARS
    yr = num2str(Y+1988-1);
    for DAY = 1:DAYS
        DAY
        n=n+1;
        %1st order upwind advect
        %bio1 = sub_advec1_diff_vel(GRD,bio1,K,U,V,ni,nj,tstep);
        bio1 = sub_diff_alt(GRD,bio1,K,ni,nj,tstep);
        biov1(:,n) = bio1(ID);
        
    end
end

% Save
save(['/Volumes/GFDL/CSV/advect_tests/Matlab_diff_K',num2str(K),'_' cname '.mat'],...
     'biov1','biov');

