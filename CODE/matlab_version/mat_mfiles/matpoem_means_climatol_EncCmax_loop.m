% POEM output at all locations
% Only last 12 months of 100 years saved to compare to SAUP for tuning

clear all
close all

encs = 10:10:100;
cmaxs = 10:10:50;
frate = 0.3; %Fish(F);

%%
for c=4:length(cmaxs)
    for e=1:length(encs)
            gam = encs(e);
            h = cmaxs(c);
            tcfn = num2str(h);
            tefn = num2str(round(gam));
            
            cfile = ['Dc_enc',tefn,'-b200_cm',tcfn,'_m-b175-k09_fcrit20',...
                '_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100'];
            fpath=['/Volumes/GFDL/CSV/Matlab_new_size/' cfile '/'];
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
            
    end
end





