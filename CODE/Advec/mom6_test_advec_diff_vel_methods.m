% Test 1st order upwind, 2nd order centered, 4th order centered 
% advection schemes from MOM5
% Input data and params needed in advection-diffusion scheme
% MOM6 grid - see if diffusion still a problem without tripolar arctic

clear all
close all

% Transports path
vpath = '/Volumes/GFDL/GCM_DATA/MOM6/Preindust/';

% Grid
load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid1D_MOM6_preindust.mat','GRD');
ID = GRD.ID;
clear GRD
load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid2D_MOM6_preindust.mat')

tmask = GRD.lmask;

%% number of water cells
NX = length(ID);

% grid size
[ni,nj] = size(tmask);
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
% cname='Global_even_dt1hr_mom6_preind_vel_surf';

%Atl-Arctic
bio(480:540,:) = 1.0e2; bio(181:241,(nj-10):nj) = 1.0e2; 
cname=['AtlArc_even_dt',tdt(2:end),'hr_mom6_preind_vel_surf'];

%seed equator
% bio(:,244:284) = 1.0e2;     
% cname='Eq_even_dt1hr_mom6_preind_vel_surf';

%seed Atl
% bio(480:540,1:500) = 1.0e2;    
% cname='Atl_even_dt1hr_mom6_preind_vel_surf';

%seed Pac
% bio(190:270,1:500) = 1.0e2;      
% cname='Pac_even_dt1hr_mom6_preind_vel_surf';

%seed Indian W
% bio(700:720,:) = 1.0e2;      
% cname='WInd_even_dt1hr_mom6_preind_vel_surf';

%seed Indian E
% bio(50:70,:) = 1.0e2;    
% cname='EInd_even_dt1hr_mom6_preind_vel_surf';

%seed Arctic
% bio(:,480:nj) = 1.0e2;    
% cname='Arc_even_dt1hr_mom6_preind_vel_surf';

%seed Antarctic
% bio(:,5:80) = 1.0e2;      
% cname='Ant_even_dt1hr_mom6_preind_vel_surf';

bio = bio .* tmask;

biov = zeros(NX,DAYS*YEARS);

%% do advec-diff
% %define diffusivity
% %horizontal diffusivity varies by two orders of magnitude from 10^2 to 10^4 m-2 s-1 
% K = 600.0;
% 
% bio1 = bio;
% bio2 = bio;
% bio4 = bio;
% biov1 = zeros(NX,DAYS*YEARS);
% biov2 = zeros(NX,DAYS*YEARS);
% biov4 = zeros(NX,DAYS*YEARS);
% B1 = zeros(ni,nj,DAYS*YEARS);
% B2 = zeros(ni,nj,DAYS*YEARS);
% B4 = zeros(ni,nj,DAYS*YEARS);
% 
% n=0;
% for Y=1:YEARS
%     % Velocities
% %     load([vpath 'MOM6_00010101_velocities_daily_interp.mat'],'D_u1','D_v1');
%     load([vpath 'MOM6_00010101_velocities_daily_interp.mat'],'D_u2','D_v2');
%     U = zeros(ni,nj);
%     V = zeros(ni,nj);
%     for DAY = 1:DAYS
%         [num2str(DAY) ',1']
%         n=n+1;
% %         u = D_u1(:,DAY);
% %         v = D_v1(:,DAY);
%         u = D_u2(:,DAY);
%         v = D_v2(:,DAY);
%         U(ID) = u; 
%         V(ID) = v;
%         %1st order upwind advect
%         bio1 = sub_mom6_advec1_diff_vel(GRD,bio1,K,U,V,ni,nj,tstep);
%         biov1(:,n) = bio1(ID);
%         B1(:,:,n) = bio1;
%         %2nd order centered advect    
% %         bio2 = sub_mom6_advec2_diff_vel(GRD,bio2,K,U,V,ni,nj,tstep);
% %         biov2(:,n) = bio2(ID);
% %         B2(:,:,n) = bio2; 
% %         %4th order centered advect
% %         bio4 = sub_mom6_advec4_diff_vel(GRD,bio4,K,U,V,ni,nj,tstep);
% %         biov4(:,n) = bio4(ID);
% %         B4(:,:,n) = bio4; 
% 
%     end
% end
% 
% % Save
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_diff_1st_2nd_4th_K',num2str(K),...
%     '_' cname '_UV2.mat'],...
%     'biov1','biov2','biov4','biov','bio2','bio4','B1','B2','B4');

%% do advec only
% K=0;
% % 1ST ORDER
% bio1 = bio;
% biov1 = zeros(NX,DAYS*YEARS);
% B1 = zeros(ni,nj,DAYS*YEARS);
% % 2ND ORDER
% bio2 = bio;
% biov2 = zeros(NX,DAYS*YEARS);
% B2 = zeros(ni,nj,DAYS*YEARS);
% % 4TH ORDER
% bio4 = bio;
% biov4 = zeros(NX,DAYS*YEARS);
% B4 = zeros(ni,nj,DAYS*YEARS);
% 
% n=0;
% for Y=1:YEARS
%     % Velocities
% %     load([vpath 'MOM6_00010101_velocities_daily_interp.mat'],'D_u1','D_v1');
%     load([vpath 'MOM6_00010101_velocities_daily_interp.mat'],'D_u2','D_v2');
%     U = zeros(ni,nj);
%     V = zeros(ni,nj);
%     for DAY = 1:DAYS
%         [num2str(DAY) ',1']
%         n=n+1;
% %         u = D_u1(:,DAY);
% %         v = D_v1(:,DAY);
%         u = D_u2(:,DAY);
%         v = D_v2(:,DAY);
%         U(ID) = u; 
%         V(ID) = v;
%         %1st order upwind advect
%         bio1 = sub_mom6_advec1_diff_vel(GRD,bio1,K,U,V,ni,nj,tstep);
%         biov1(:,n) = bio1(ID);
%         B1(:,:,n) = bio1;
%         %2nd order centered advect    
%         bio2 = sub_mom6_advec2_diff_vel(GRD,bio2,K,U,V,ni,nj,tstep);
%         biov2(:,n) = bio2(ID);
%         B2(:,:,n) = bio2; 
%         %4th order centered advect
%         bio4 = sub_mom6_advec4_diff_vel(GRD,bio4,K,U,V,ni,nj,tstep);
%         biov4(:,n) = bio4(ID);
%         B4(:,:,n) = bio4; 
% 
%     end
% end
% 
% % Save
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_1st_2nd_4th_' cname '_UV2.mat'],...
%     'biov1','biov2','biov4','biov','bio2','bio4','B1','B2','B4');

%% do diff only
%Global
% bio = 100*rand(ni,nj);   
% cname='Global_rand_dt1hr_mom6_preind_vel_surf';
%Gradient highest at eq
% g1 = 1:100;
% g2 = 100:-1:1;
% bio = zeros(ni,nj);
% bio(:,1:100) = repmat(g1,ni,1);
% bio(:,101:200) = repmat(g2,ni,1);
% cname='Eq_grad_dt1hr_mom6_preind_vel_surf';
%Gradient highest at poles
% g1 = 1:100;
% g2 = 100:-1:1;
% bio = zeros(ni,nj);
% bio(:,1:100) = repmat(g2,ni,1);
% bio(:,477:576) = repmat(g1,ni,1);
% cname='Pole_grad_dt1hr_mom6_preind_vel_surf';

bio = bio .* GRD.lmask;

bio1 = bio;
bio2 = bio;
bio4 = bio;
biov1 = zeros(NX,DAYS*YEARS);
biov2 = zeros(NX,DAYS*YEARS);
biov4 = zeros(NX,DAYS*YEARS);
B1 = zeros(ni,nj,DAYS*YEARS);
B2 = zeros(ni,nj,DAYS*YEARS);
B4 = zeros(ni,nj,DAYS*YEARS);

% define diffusivity
K = 1;
U = zeros(ni,nj);
V = zeros(ni,nj);
n=0;
for Y=1:YEARS
    for DAY = 1:DAYS
        [num2str(DAY) ',1']
        n=n+1;
        % all versions have same diffusion
        %1st order upwind advect
        bio1 = sub_mom6_advec1_diff_vel(GRD,bio1,K,U,V,ni,nj,tstep);
        biov1(:,n) = bio1(ID);
        B1(:,:,n) = bio1;
%         %2nd order centered advect    
%         bio2 = sub_mom6_advec2_diff_vel(GRD,bio2,K,U,V,ni,nj,tstep);
%         biov2(:,n) = bio2(ID);
%         B2(:,:,n) = bio2; 
%         %4th order centered advect
%         bio4 = sub_mom6_advec4_diff_vel(GRD,bio4,K,U,V,ni,nj,tstep);
%         biov4(:,n) = bio4(ID);
%         B4(:,:,n) = bio4; 
    end
end
% Save
save(['/Volumes/GFDL/CSV/advect_tests/Matlab_diff_K',num2str(K),...
    '_' cname '.mat'],...
     'biov1','bio1','B1');

