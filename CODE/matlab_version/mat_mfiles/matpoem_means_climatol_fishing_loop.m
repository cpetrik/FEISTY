% POEM output at all locations

clear all
close all

fqs = 0.5:0.5:3.5;
pqs = 0.25:0.25:1.5;
dqs = 1:7;
frate = 0.1; %Fish(F);
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE10_CC100_lgRE00100_mdRE00100';
fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];

%%
% for fq=1:length(fqs)
%     for pq=4;%1:length(pqs)
%         for dq=7;%1:length(dqs)
%             MFsel = fqs(fq);
%             LPsel = pqs(pq);
%             LDsel = dqs(dq);
%             tF = num2str(1000+int64(100*frate*MFsel));
%             tP = num2str(1000+int64(100*frate*LPsel));
%             tD = num2str(1000+int64(100*frate*LDsel));
            
%             sfile = ['/Volumes/GFDL/CSV/Matlab_new_size/',cfile,...
%                 '/Clim_means_fish_qF',tF(2:end),'_qP',tP(2:end),'_qD',tD(2:end),'.mat'];
%             load(sfile);
%             harv = ['fish_qF',tF(2:end),'_qP',tP(2:end),'_qD',tD(2:end)];
%             
            sfile = ['/Volumes/GFDL/CSV/Matlab_new_size/',cfile,...
                '/Clim_means_All_fish03.mat'];
            load(sfile);
            harv = 'All_fish03';
            
            % Last year
            [id,nt] = size(Spinup_Bent.bio);
            time=1:nt;
            lyr=time((end-12+1):end);
            sp_mean=mean(Spinup_Sml_p.bio(:,lyr),2);
            sf_mean=mean(Spinup_Sml_f.bio(:,lyr),2);
            sd_mean=mean(Spinup_Sml_d.bio(:,lyr),2);
            mp_mean=mean(Spinup_Med_p.bio(:,lyr),2);
            mf_mean=mean(Spinup_Med_f.bio(:,lyr),2);
            md_mean=mean(Spinup_Med_d.bio(:,lyr),2);
            lp_mean=mean(Spinup_Lrg_p.bio(:,lyr),2);
            ld_mean=mean(Spinup_Lrg_d.bio(:,lyr),2);
            b_mean=mean(Spinup_Bent.bio(:,lyr),2);
            
            mf_my=mean(Spinup_Med_f.yield(:,lyr),2);
            mp_my=mean(Spinup_Med_p.yield(:,lyr),2);
            md_my=mean(Spinup_Med_d.yield(:,lyr),2);
            lp_my=mean(Spinup_Lrg_p.yield(:,lyr),2);
            ld_my=mean(Spinup_Lrg_d.yield(:,lyr),2);
            
            %
            save(sfile,...
                'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
                'lp_mean','ld_mean','time','lyr',...
                'mf_my','mp_my','md_my','lp_my','ld_my',...
                '-append');
            
%         end
%     end
% end





