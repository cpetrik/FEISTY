% Make 3 or 4-D matrices of Climatol fishing diff F results

clear all
close all

fpath = '/Volumes/FEISTY/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/Climatology/';

%%
load([fpath 'Climatol_fish_diffF_tot_catch1.mat']);
load('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF.mat');

frates = [0:0.1:2]';

fstr = num2str(1000+int64(100*frates));
fstr = fstr(:,2:end);
fstr2 = cellstr(fstr);

%%
Fy = nan*ones(length(frates),length(frates),length(frates));
Py = Fy;
Dy = Fy;

for i=1:length(frates)
    for j=1:length(frates)
        for k=1:length(frates)
            fF = fstr2{i};
            fP = fstr2{j};
            fD = fstr2{k};
            
            rname = ['Climatol_fish_F',fF,'_P',fP,'_D',fD,'.mat'];
            if strcmp(rname,'Climatol_fish_F000_P000_D000.mat')
                rname = 'Climatol_pristine';
            end
            
            lid = strcmp(rname,Fnames);
            fid = find(lid==1);
            if (~isempty(fid))
                Fy(i,j,k) = mf_gtc(fid);
                Py(i,j,k) = mp_gtc(fid) + lp_gtc(fid);
                Dy(i,j,k) = md_gtc(fid) + ld_gtc(fid);
            end
            clear fid lid
            
        end
    end
end

%%
save([fpath 'Climatol_fish_diffF_tot_catch1.mat'],'Fy','Py','Dy',...
    '-append');

%% test viz
[x,y] = meshgrid(frates,frates);

figure(1)
subplot(2,2,1)
pcolor(x,y,squeeze(Fy(:,:,1)))
colorbar
title('F yield F & P')

subplot(2,2,2)
pcolor(x,y,squeeze(Fy(:,1,:)))
colorbar
title('F yield F & D')

subplot(2,2,3)
pcolor(x,y,squeeze(Fy(1,:,:)))
colorbar
title('F yield P & D')


figure(2)
subplot(2,2,1)
pcolor(x,y,squeeze(Py(:,:,1)))
colorbar
title('P yield F & P')

subplot(2,2,2)
pcolor(x,y,squeeze(Py(:,1,:)))
colorbar
title('P yield F & D')

subplot(2,2,3)
pcolor(x,y,squeeze(Py(1,:,:)))
colorbar
title('P yield P & D')


figure(3)
subplot(2,2,1)
pcolor(x,y,squeeze(Dy(:,:,1)))
colorbar
title('D yield F & P')

subplot(2,2,2)
pcolor(x,y,squeeze(Dy(:,1,:)))
colorbar
title('D yield F & D')

subplot(2,2,3)
pcolor(x,y,squeeze(Dy(1,:,:)))
colorbar
title('D yield P & D')

