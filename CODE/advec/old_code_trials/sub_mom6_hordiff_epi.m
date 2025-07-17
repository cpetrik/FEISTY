%> Main routine for lateral (along surface or neutral) diffusion of tracers
%module MOM_tracer_hor_diff

% This file is part of MOM6. See LICENSE.md for the license.
%use MOM_neutral_diffusion.F90

%> Compute along-coordinate diffusion of all tracers
% using the diffusivity in CS%KhTr, or using space-dependent diffusivity.
% Multiple iterations are used (if necessary) so that there is no limit
% on the acceptable time increment.

% ===================================================================================
% subroutine tracer_hordiff(h, dt, MEKE, VarMix, G, GV, CS, Reg, tv, do_online_flag, read_khdt_x, read_khdt_y)
%   G       %< Grid type
%   h       %< Layer thickness [H ~> m or kg m-2]
%   dt      %< time step [s]
%   CS      %< module control structure
%   Ihdxdy, &     % The inverse of the volume or mass of fluid in a layer in a
%                   % grid cell [H-1 m-2 ~> m-3 or kg-1].
%     Kh_h, &       % The tracer diffusivity averaged to tracer points [m2 s-1].
%     CFL, &        % A diffusive CFL number for each cell [nondim].
%     dTr           % The change in a tracer's concentration, in units of concentration [Conc].
%     khdt_x, &     % The value of Khtr*dt times the open face width divided by
%                   % the distance between adjacent tracer points [m2].
%     Coef_x, &     % The coefficients relating zonal tracer differences
%                   % to time-integrated fluxes [H m2 ~> m3 or kg].
%     Kh_u          % Tracer mixing coefficient at u-points [m2 s-1].
%     khdt_y, &     % The value of Khtr*dt times the open face width divided by
%                   % the distance between adjacent tracer points [m2].
%     Coef_y, &     % The coefficients relating meridional tracer differences
%                   % to time-integrated fluxes [H m2 ~> m3 or kg].
%     Kh_v          % Tracer mixing coefficient at u-points [m2 s-1].
%
%   khdt_max % The local limiting value of khdt_x or khdt_y [m2].
%   max_CFL % The global maximum of the diffusive CFL number.
%   I_numitts  % The inverse of the number of iterations, num_itts.
%   scale      % The fraction of khdt_x or khdt_y that is applied in this
%                      % layer for this iteration [nondim].
%   Idt        % The inverse of the time step [s-1].
%   h_neglect  % A thickness that is so small it is usually lost
%                      % in roundoff and can be neglected [H ~> m or kg m-2].

clear all
close all

% Transports path
vpath = '/Volumes/GFDL/GCM_DATA/MOM6/Preindust/';

% Grid
load('/Volumes/GFDL/GCM_DATA/MOM6/Data_grid2D_diff_MOM6_preindust.mat')

%% define time
tstep = 1; %time step in hours
dt = 60*60*tstep; %time step in sec

%% define a patch to advect
[ni,nj] = size(G.area);
bio = zeros(ni,nj);
%File name to save

%Global
% bio = 100*ones(ni,nj);   %Global
% cname='Global_even_dt1hr_mom6_preind_vel_surf';

%Atl-Arctic
bio(480:540,:) = 1.0e2; bio(181:241,(nj-10):nj) = 1.0e2; 
cname=['AtlArc_even_dt',num2str(dt),'s_mom6_preind'];

Tr = bio;

%%
is = G.isc ;
ie = G.iec ;
js = G.jsc ;
je = G.jec ;
%nz = GV.ke;

Idt = 1.0/dt;
h_neglect = eps; % A thickness that is so small it is usually lost in roundoff and can be neglected [H ~> m or kg m-2].

CS.max_diff_CFL = 1.0;
CS.use_neutral_diffusion = 0;
CS.Diffuse_ML_interior = 1;

% Use a simple constant diffusivity.
K = 600;
h = G.mask; %layer thickness -> assume 1 m because conc in g/m2 ?
G.IdxCu = 1./G.dxCu;
G.IdyCv = 1./G.dyCv;
khdt_x = dt*K*(G.dyCu .* G.IdxCu);
khdt_y = dt*K*(G.dxCv .* G.IdyCv);

%%
if (CS.max_diff_CFL > 0.0)
    for j=js:je
        %for i=is-1:ie
        for i=is:ie
            if (i==ie)
                khdt_max = 0.125*CS.max_diff_CFL * min(G.area(i,j), G.area(is,j));
            else
                khdt_max = 0.125*CS.max_diff_CFL * min(G.area(i,j), G.area(i+1,j));
            end
            khdt_x(i,j) = min(khdt_x(i,j), khdt_max);
        end
    end
    
    for j=js-1:je
        for i=is:ie
            if (j==je)
                khdt_max = 0.125*CS.max_diff_CFL * min(G.area(i,j), G.area(ni-i+1,j));
            else
                khdt_max = 0.125*CS.max_diff_CFL * min(G.area(i,j), G.area(i,j+1));
            end
            khdt_y(i,j) = min(khdt_y(i,j), khdt_max);
        end
    end
end

if (CS.max_diff_CFL > 0.0)
    %num_itts = max(1, ceiling(CS.max_diff_CFL - 4.0*epsilon(CS.max_diff_CFL)));
    num_itts = max(1, ceil(CS.max_diff_CFL - 4.0*eps*CS.max_diff_CFL));
    I_numitts = 1.0 / (real(num_itts));
else
    num_itts = 1;
    I_numitts = 1.0;
end
%%
Coef_x = nan(ni,nj);
Coef_y = nan(ni,nj);
Ihdxdy = nan(ni,nj);
dTr    = nan(ni,nj);

if (CS.use_neutral_diffusion)
    %call neutral_diffusion_calc_coeffs(G, GV, h, tv%T, tv%S, CS.neutral_diffusion_CSp)
    Coef_y = I_numitts .* khdt_y;
    Coef_x = I_numitts .* khdt_x;
    % Coef_x The coefficients relating zonal tracer differences to time-integrated fluxes [H m2 ~> m3 or kg].
    % Coef_y The coefficients relating meridional tracer differences to time-integrated fluxes [H m2 ~> m3 or kg].
    for itt=1:num_itts
        Tr = neutral_diffusion(G,GV,h,Coef_x,Coef_y,I_numitts*dt,Tr);
    end % itt
    
else  % not using neutral diffusion, but instead along-surface diffusion
    
    for itt=1:num_itts
        scale = I_numitts;
        
        for j=js-1:je
            for i=is:ie
                if (j==je)
                    Coef_y(i,j) = ((scale * khdt_y(i,j))*2.0*(h(i,j)*h(ni-i+1,j))) / ...
                    (h(i,j)+h(ni-i+1,j)+h_neglect);
                else
                    Coef_y(i,j) = ((scale * khdt_y(i,j))*2.0*(h(i,j)*h(i,j+1))) / ...
                    (h(i,j)+h(i,j+1)+h_neglect);
                end
            end
        end
        
        for j=js:je
            %for i=is-1:ie
            for i=is:ie
                if (i==ie)
                    Coef_x(i,j) = ((scale * khdt_x(i,j))*2.0*(h(i,j)*h(is,j))) / ...
                    (h(i,j)+h(is,j)+h_neglect);
                else
                    Coef_x(i,j) = ((scale * khdt_x(i,j))*2.0*(h(i,j)*h(i+1,j))) / ...
                    (h(i,j)+h(i+1,j)+h_neglect);
                end
            end
            for i=is:ie
                Ihdxdy(i,j) = G.Iarea(i,j) / (h(i,j)+h_neglect);
            end
        end
        
        for j=js:je
            for i=is:ie
                if (j==je)
                    if (i==ie)
                        dTr(i,j) = Ihdxdy(i,j) * ...
                        ((Coef_x(i-1,j) * (Tr(i-1,j) - Tr(i,j)) - ...
                        Coef_x(i,j) * (Tr(i,j) - Tr(is,j))) + ...
                        (Coef_y(i,j-1) * (Tr(i,j-1) - Tr(i,j)) - ...
                        Coef_y(i,j) * (Tr(i,j) - Tr(ni-i+1,j))));
                    elseif (i==is)
                        dTr(i,j) = Ihdxdy(i,j) * ...
                        ((Coef_x(ie,j) * (Tr(ie,j) - Tr(i,j)) - ...
                        Coef_x(i,j) * (Tr(i,j) - Tr(i+1,j))) + ...
                        (Coef_y(i,j-1) * (Tr(i,j-1) - Tr(i,j)) - ...
                        Coef_y(i,j) * (Tr(i,j) - Tr(ni-i+1,j))));
                    else
                        dTr(i,j) = Ihdxdy(i,j) * ...
                        ((Coef_x(i-1,j) * (Tr(i-1,j) - Tr(i,j)) - ...
                        Coef_x(i,j) * (Tr(i,j) - Tr(i+1,j))) + ...
                        (Coef_y(i,j-1) * (Tr(i,j-1) - Tr(i,j)) - ...
                        Coef_y(i,j) * (Tr(i,j) - Tr(ni-i+1,j))));
                    end
                else
                    if (i==ie)
                        dTr(i,j) = Ihdxdy(i,j) * ...
                            ((Coef_x(i-1,j) * (Tr(i-1,j) - Tr(i,j)) - ...
                            Coef_x(i,j) * (Tr(i,j) - Tr(is,j))) + ...
                            (Coef_y(i,j-1) * (Tr(i,j-1) - Tr(i,j)) - ...
                            Coef_y(i,j) * (Tr(i,j) - Tr(i,j+1))));
                    elseif (i==is)
                        dTr(i,j) = Ihdxdy(i,j) * ...
                            ((Coef_x(ie,j) * (Tr(ie,j) - Tr(i,j)) - ...
                            Coef_x(i,j) * (Tr(i,j) - Tr(i+1,j))) + ...
                            (Coef_y(i,j-1) * (Tr(i,j-1) - Tr(i,j)) - ...
                            Coef_y(i,j) * (Tr(i,j) - Tr(i,j+1))));
                    else
                        dTr(i,j) = Ihdxdy(i,j) * ...
                            ((Coef_x(i-1,j) * (Tr(i-1,j) - Tr(i,j)) - ...
                            Coef_x(i,j) * (Tr(i,j) - Tr(i+1,j))) + ...
                            (Coef_y(i,j-1) * (Tr(i,j-1) - Tr(i,j)) - ...
                            Coef_y(i,j) * (Tr(i,j) - Tr(i,j+1))));
                    end
                end
            end
        end
        
%         if (~isnan(Reg.df_x))
%             for j=js,je ;
%                 for i=G.IscB:G.IecB
%                     Reg.df_x(i,j) = Reg.df_x(i,j) + Coef_x(i,j) * ...
%                         (Tr(i,j) - Tr(i+1,j))*Idt;
%                 end
%             end
%         end
%         
%         if (associated(Reg.df_y))
%             for j=G.JscB,G.JecB ;
%                 for i=is:ie
%                     Reg.df_y(i,j) = Reg.df_y(i,j) + Coef_y(i,j) * ...
%                         (Tr(i,j) - Tr(i,j+1))*Idt;
%                 end
%             end
%         end
%         
%         if (associated(Reg.df2d_x))
%             for j=js,je ;
%                 for i=G.IscB,G.IecB
%                     Reg.df2d_x(i,j) = Reg.df2d_x(i,j) + Coef_x(i,j) * ...
%                         (Tr(i,j) - Tr(i+1,j))*Idt;
%                 end
%             end
%         end
%         
%         if (associated(Reg.df2d_y))
%             for j=G.JscB:G.JecB ;
%                 for i=is:ie
%                     Reg.df2d_y(i,j) = Reg.df2d_y(i,j) + Coef_y(i,j) * ...
%                         (Tr(i,j) - Tr(i,j+1))*Idt;
%                 end
%             end
%         end
        
        for j=js:je
            for i=is:ie
                Tr(i,j) = Tr(i,j) + dTr(i,j);
            end
        end
    end % End of nt loop.
    
end   % endif for CS.use_neutral_diffusion

%%
if (CS.Diffuse_ML_interior)
    Tr = sub_mom6_epipycnal_diff(h, dt, Tr, K, G, num_itts);
end


%end subroutine tracer_hordiff ==========================================

%
%    This module contains subroutines that handle horizontal
%  diffusion (i.e., isoneutral or along layer) of tracers.
%
%    Each of the tracers are subject to Fickian along-coordinate
%  diffusion if Khtr is defined and positive.  The tracer diffusion
%  can use a suitable number of iterations to guarantee stability
%  with an arbitrarily large time step.

%end module MOM_tracer_hor_diff


