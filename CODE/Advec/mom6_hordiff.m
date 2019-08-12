%> Main routine for lateral (along surface or neutral) diffusion of tracers
%module MOM_tracer_hor_diff

% This file is part of MOM6. See LICENSE.md for the license.
%use MOM_neutral_diffusion.F90

%% > Compute along-coordinate diffusion of all tracers
% using the diffusivity in CS%KhTr, or using space-dependent diffusivity.
% Multiple iterations are used (if necessary) so that there is no limit
% on the acceptable time increment.

% ===================================================================================
% subroutine tracer_hordiff(h, dt, MEKE, VarMix, G, GV, CS, Reg, tv, do_online_flag, read_khdt_x, read_khdt_y)
%   type(ocean_grid_type),      intent(inout) :: G       %< Grid type
%   real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
%                               intent(in)    :: h       %< Layer thickness [H ~> m or kg m-2]
%   real,                       intent(in)    :: dt      %< time step [s]
%   type(MEKE_type),            pointer       :: MEKE    %< MEKE type
%   type(VarMix_CS),            pointer       :: VarMix  %< Variable mixing type
%   type(verticalGrid_type),    intent(in)    :: GV      %< ocean vertical grid structure
%   type(tracer_hor_diff_CS),   pointer       :: CS      %< module control structure
%   type(tracer_registry_type), pointer       :: Reg     %< registered tracers
%   type(thermo_var_ptrs),      intent(in)    :: tv      %< A structure containing pointers to any available
%                                                        %% thermodynamic fields, including potential temp and
%                                                        %% salinity or mixed layer density. Absent fields have
%                                                        %% NULL ptrs, and these may (probably will) point to
%                                                        %% some of the same arrays as Tr does.  tv is required
%                                                        %% for epipycnal mixing between mixed layer and the interior.
%   % Optional inputs for offline tracer transport
%   logical,          optional, intent(in)    :: do_online_flag %< If present and true, do online
%                                                        %% tracer transport with stored velcities.
%   real, dimension(SZIB_(G),SZJ_(G)), &
%                     optional, intent(in)    :: read_khdt_x %< If present, these are the zonal
%                                                        %% diffusivities from previous run.
%   real, dimension(SZI_(G),SZJB_(G)), &
%                     optional, intent(in)    :: read_khdt_y %< If present, these are the meridional
%                                                        %% diffusivities from previous run.
% 
% 
%   real, dimension(SZI_(G),SZJ_(G)) :: &
%     Ihdxdy, &     % The inverse of the volume or mass of fluid in a layer in a
%                   % grid cell [H-1 m-2 ~> m-3 or kg-1].
%     Kh_h, &       % The tracer diffusivity averaged to tracer points [m2 s-1].
%     CFL, &        % A diffusive CFL number for each cell [nondim].
%     dTr           % The change in a tracer's concentration, in units of concentration [Conc].
% 
%   real, dimension(SZIB_(G),SZJ_(G)) :: &
%     khdt_x, &     % The value of Khtr*dt times the open face width divided by
%                   % the distance between adjacent tracer points [m2].
%     Coef_x, &     % The coefficients relating zonal tracer differences
%                   % to time-integrated fluxes [H m2 ~> m3 or kg].
%     Kh_u          % Tracer mixing coefficient at u-points [m2 s-1].
%   real, dimension(SZI_(G),SZJB_(G)) :: &
%     khdt_y, &     % The value of Khtr*dt times the open face width divided by
%                   % the distance between adjacent tracer points [m2].
%     Coef_y, &     % The coefficients relating meridional tracer differences
%                   % to time-integrated fluxes [H m2 ~> m3 or kg].
%     Kh_v          % Tracer mixing coefficient at u-points [m2 s-1].
% 
%   real :: khdt_max % The local limiting value of khdt_x or khdt_y [m2].
%   real :: max_CFL % The global maximum of the diffusive CFL number.
%   logical :: use_VarMix, Resoln_scaled, do_online, use_Eady
%   integer :: S_idx, T_idx % Indices for temperature and salinity if needed
%   integer :: i, j, k, m, is, ie, js, je, nz, ntr, itt, num_itts
%   real :: I_numitts  % The inverse of the number of iterations, num_itts.
%   real :: scale      % The fraction of khdt_x or khdt_y that is applied in this
%                      % layer for this iteration [nondim].
%   real :: Idt        % The inverse of the time step [s-1].
%   real :: h_neglect  % A thickness that is so small it is usually lost
%                      % in roundoff and can be neglected [H ~> m or kg m-2].
%   real :: Kh_loc     % The local value of Kh [m2 s-1].
%   real :: Res_Fn     % The local value of the resolution function [nondim].
%   real :: Rd_dx      % The local value of deformation radius over grid-spacing [nondim].
%   real :: normalize  % normalization used for diagnostic Kh_h; diffusivity averaged to h-points.

  is = GRD.isc ; ie = GRD.iec ; js = GRD.jsc ; je = GRD.jec ; nz = GV.ke;

  Idt = 1.0/dt
  %h_neglect = GV.H_subroundoff

  
  % Use a simple constant diffusivity.
      if (CS.id_KhTr_u > 0) then
        %$OMP parallel do default(shared)
        do j=js,je ; do I=is-1,ie
          Kh_u(I,j) = CS.KhTr
          khdt_x(I,j) = dt*(CS.KhTr*(GRD.dy_Cu(I,j)*GRD.IdxCu(I,j)))
        enddo ; enddo
      else
        %$OMP parallel do default(shared)
        do j=js,je ; do I=is-1,ie
          khdt_x(I,j) = dt*(CS.KhTr*(GRD.dy_Cu(I,j)*GRD.IdxCu(I,j)))
        enddo ; enddo
      end
      if (CS.id_KhTr_v > 0) then
        %$OMP parallel do default(shared)
        do J=js-1,je ;  do i=is,ie
          Kh_v(i,J) = CS.KhTr
          khdt_y(i,J) = dt*(CS.KhTr*(GRD.dx_Cv(i,J)*GRD.IdyCv(i,J)))
        enddo ; enddo
      else
        %$OMP parallel do default(shared)
        do J=js-1,je ;  do i=is,ie
          khdt_y(i,J) = dt*(CS.KhTr*(GRD.dx_Cv(i,J)*GRD.IdyCv(i,J)))
        enddo ; enddo
      end
    

    if (CS.max_diff_CFL > 0.0) then
      if ((CS.id_KhTr_u > 0) .or. (CS.id_KhTr_h > 0)) then
        %$OMP parallel do default(shared) private(khdt_max)
        do j=js,je ; do I=is-1,ie
          khdt_max = 0.125*CS.max_diff_CFL * min(GRD.areaT(i,j), GRD.areaT(i+1,j))
          if (khdt_x(I,j) > khdt_max) then
            khdt_x(I,j) = khdt_max
            if (dt*(GRD.dy_Cu(I,j)*GRD.IdxCu(I,j)) > 0.0) &
              Kh_u(I,j) = khdt_x(I,j) / (dt*(GRD.dy_Cu(I,j)*GRD.IdxCu(I,j)))
            end
        enddo ; enddo
      else
        %$OMP parallel do default(shared) private(khdt_max)
        do j=js,je ; do I=is-1,ie
          khdt_max = 0.125*CS.max_diff_CFL * min(GRD.areaT(i,j), GRD.areaT(i+1,j))
          khdt_x(I,j) = min(khdt_x(I,j), khdt_max)
        enddo ; enddo
      end
      if ((CS.id_KhTr_v > 0) .or. (CS.id_KhTr_h > 0)) then
        %$OMP parallel do default(shared) private(khdt_max)
        do J=js-1,je ; do i=is,ie
          khdt_max = 0.125*CS.max_diff_CFL * min(GRD.areaT(i,j), GRD.areaT(i,j+1))
          if (khdt_y(i,J) > khdt_max) then
            khdt_y(i,J) = khdt_max
            if (dt*(GRD.dx_Cv(i,J)*GRD.IdyCv(i,J)) > 0.0) &
              Kh_v(i,J) = khdt_y(i,J) / (dt*(GRD.dx_Cv(i,J)*GRD.IdyCv(i,J)))
          end
        enddo ; enddo
      else
        %$OMP parallel do default(shared) private(khdt_max)
        do J=js-1,je ; do i=is,ie
          khdt_max = 0.125*CS.max_diff_CFL * min(GRD.areaT(i,j), GRD.areaT(i,j+1))
          khdt_y(i,J) = min(khdt_y(i,J), khdt_max)
        enddo ; enddo
      end
    end

  if (CS.check_diffusive_CFL) then
    max_CFL = 0.0
    do j=js,je ; do i=is,ie
      CFL(i,j) = 2.0*((khdt_x(I-1,j) + khdt_x(I,j)) + &
                      (khdt_y(i,J-1) + khdt_y(i,J))) * GRD.IareaT(i,j)
      if (max_CFL < CFL(i,j)) max_CFL = CFL(i,j)
    enddo ; enddo
    num_itts = max(1, ceiling(max_CFL - 4.0*EPSILON(max_CFL)))
    I_numitts = 1.0 / (real(num_itts))
    if (CS.id_CFL > 0) call post_data(CS.id_CFL, CFL, CS.diag, mask=GRD.mask2dT)
  elseif (CS.max_diff_CFL > 0.0) then
    num_itts = max(1, ceiling(CS.max_diff_CFL - 4.0*EPSILON(CS.max_diff_CFL)))
    I_numitts = 1.0 / (real(num_itts))
  else
    num_itts = 1 ; I_numitts = 1.0
    end
    

  if (CS.use_neutral_diffusion) then

    %call neutral_diffusion_calc_coeffs(G, GV, h, tv%T, tv%S, CS.neutral_diffusion_CSp)
    do J=js-1,je ; do i=is,ie
      Coef_y(i,J) = I_numitts * khdt_y(i,J)
    enddo ; enddo
    do j=js,je
      do I=is-1,ie
        Coef_x(I,j) = I_numitts * khdt_x(I,j)
      enddo
    enddo

%     Coef_x, &     % The coefficients relating zonal tracer differences
%                   % to time-integrated fluxes [H m2 ~> m3 or kg].
%     Coef_y, &     % The coefficients relating meridional tracer differences
%                   % to time-integrated fluxes [H m2 ~> m3 or kg].
    do itt=1,num_itts
      call neutral_diffusion(G, GV,  h, Coef_x, Coef_y, I_numitts*dt, Reg, CS.neutral_diffusion_CSp)
    enddo % itt

  end   % end for CS.use_neutral_diffusion
  

%end subroutine tracer_hordiff ==========================================

% ====================================================================================
%> Calculate remapping factors for u/v columns used to map adjoining columns to
%% a shared coordinate space.
%%subroutine neutral_diffusion_calc_coeffs(G, GV, h, T, S, CS)
%   type(ocean_grid_type),                    intent(in) :: G   %< Ocean grid structure
%   type(verticalGrid_type),                  intent(in) :: GV  %< ocean vertical grid structure
%   real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h   %< Layer thickness [H ~> m or kg m-2]
%   real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: T   %< Potential temperature [degC]
%   real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: S   %< Salinity [ppt]
%   type(neutral_diffusion_CS),               pointer    :: CS  %< Neutral diffusion control structure
% 
%   % Local variables
%   integer :: i, j, k
%   % Variables used for reconstructions
%   real, dimension(SZK_(G),2) :: ppoly_r_S            % Reconstruction slopes
%   real, dimension(SZI_(G), SZJ_(G)) :: hEff_sum
%   integer :: iMethod
%   real, dimension(SZI_(G)) :: ref_pres % Reference pressure used to calculate alpha/beta
%   real :: h_neglect, h_neglect_edge
%   real :: pa_to_H
% 
%   pa_to_H = 1. / GV.H_to_pa
% 
%   %### Try replacing both of these with GV.H_subroundoff
%   h_neglect = GV.kg_m2_to_H*1.0e-30 ; 
%   h_neglect_edge = GV.kg_m2_to_H*1.0e-10
%   
%   if (CS.continuous_reconstruction) then
%     CS.dRdT(:,:,:) = 0.
%     CS.dRdS(:,:,:) = 0.
%   else
%     CS.T_i(:,:,:,:) = 0.
%     CS.S_i(:,:,:,:) = 0.
%     CS.dRdT_i(:,:,:,:) = 0.
%     CS.dRdS_i(:,:,:,:) = 0.
%     CS.ns(:,:) = 0.
%     CS.stable_cell(:,:,:) = .true.
%   end
% 
%   % Calculate pressure at interfaces and layer averaged alpha/beta
%   CS.Pint(:,:,1) = 0.
%   do k=1,GRD.ke ; do j=GRD.jsc-1, GRD.jec+1 ; do i=GRD.isc-1,GRD.iec+1
%       CS.Pint(i,j,k+1) = CS.Pint(i,j,k) + h(i,j,k)*GV.H_to_Pa
%   enddo ; enddo ; enddo
% 
%   do j = GRD.jsc-1, GRD.jec+1
%     % Interpolate state to interface
%     do i = GRD.isc-1, GRD.iec+1
%       if (CS.continuous_reconstruction) then
%         call interface_scalar(GRD.ke, h(i,j,:), T(i,j,:), CS.Tint(i,j,:), 2, h_neglect)
%         call interface_scalar(GRD.ke, h(i,j,:), S(i,j,:), CS.Sint(i,j,:), 2, h_neglect)
%       else
%         call build_reconstructions_1d( CS.remap_CS, GRD.ke, h(i,j,:), T(i,j,:), CS.ppoly_coeffs_T(i,j,:,:), &
%                                        CS.T_i(i,j,:,:), ppoly_r_S, iMethod, h_neglect, h_neglect_edge )
%         call build_reconstructions_1d( CS.remap_CS, GRD.ke, h(i,j,:), S(i,j,:), CS.ppoly_coeffs_S(i,j,:,:), &
%                                        CS.S_i(i,j,:,:), ppoly_r_S, iMethod, h_neglect, h_neglect_edge )
%       end
%     enddo
% 
%     % Continuous reconstruction
%     if (CS.continuous_reconstruction) then
%       do k = 1, GRD.ke+1
%         if (CS.ref_pres<0) ref_pres(:) = CS.Pint(:,j,k)
%         call calculate_density_derivs(CS.Tint(:,j,k), CS.Sint(:,j,k), ref_pres, &
%                                       CS.dRdT(:,j,k), CS.dRdS(:,j,k), GRD.isc-1, GRD.iec-GRD.isc+3, CS.EOS)
%       enddo
%     else % Discontinuous reconstruction
%       do k = 1, GRD.ke
%         if (CS.ref_pres<0) ref_pres(:) = CS.Pint(:,j,k)
%         call calculate_density_derivs(CS.T_i(:,j,k,1), CS.S_i(:,j,k,1), ref_pres, &
%                                       CS.dRdT_i(:,j,k,1), CS.dRdS_i(:,j,k,1), GRD.isc-1, GRD.iec-GRD.isc+3, CS.EOS)
%         if (CS.ref_pres<0) then
%           ref_pres(:) = CS.Pint(:,j,k+1)
%         end
%         call calculate_density_derivs(CS.T_i(:,j,k,2), CS.S_i(:,j,k,2), ref_pres, &
%                                       CS.dRdT_i(:,j,k,2), CS.dRdS_i(:,j,k,2), GRD.isc-1, GRD.iec-GRD.isc+3, CS.EOS)
%       enddo
%     end
%   enddo
% 
%   if (.not. CS.continuous_reconstruction) then
%     do j = GRD.jsc-1, GRD.jec+1 ; do i = GRD.isc-1, GRD.iec+1
%       call mark_unstable_cells( GRD.ke, CS.dRdT_i(i,j,:,:), CS.dRdS_i(i,j,:,:), CS.T_i(i,j,:,:), CS.S_i(i,j,:,:), &
%                                 CS.stable_cell(i,j,:), CS.ns(i,j) )
%     enddo ; enddo
%   end
% 
%   CS.uhEff(:,:,:) = 0.
%   CS.vhEff(:,:,:) = 0.
%   CS.uPoL(:,:,:) = 0.
%   CS.vPoL(:,:,:) = 0.
%   CS.uPoR(:,:,:) = 0.
%   CS.vPoR(:,:,:) = 0.
%   CS.uKoL(:,:,:) = 1
%   CS.vKoL(:,:,:) = 1
%   CS.uKoR(:,:,:) = 1
%   CS.vKoR(:,:,:) = 1
% 
%   % Neutral surface factors at U points
%   do j = GRD.jsc, GRD.jec ; do I = GRD.isc-1, GRD.iec
%     if (GRD.mask2dCu(I,j) > 0.) then
%       if (CS.continuous_reconstruction) then
%         call find_neutral_surface_positions_continuous(GRD.ke,                                    &
%                 CS.Pint(i,j,:), CS.Tint(i,j,:), CS.Sint(i,j,:), CS.dRdT(i,j,:), CS.dRdS(i,j,:),            &
%                 CS.Pint(i+1,j,:), CS.Tint(i+1,j,:), CS.Sint(i+1,j,:), CS.dRdT(i+1,j,:), CS.dRdS(i+1,j,:),  &
%                 CS.uPoL(I,j,:), CS.uPoR(I,j,:), CS.uKoL(I,j,:), CS.uKoR(I,j,:), CS.uhEff(I,j,:) )
%       else
%         call find_neutral_surface_positions_discontinuous(CS, GRD.ke, CS.ns(i,j)+CS.ns(i+1,j), &
%             CS.Pint(i,j,:), h(i,j,:), CS.T_i(i,j,:,:), CS.S_i(i,j,:,:),                      &
%             CS.dRdT_i(i,j,:,:), CS.dRdS_i(i,j,:,:), CS.stable_cell(i,j,:),                   &
%             CS.Pint(i+1,j,:), h(i+1,j,:), CS.T_i(i+1,j,:,:), CS.S_i(i+1,j,:,:),              &
%             CS.dRdT_i(i+1,j,:,:), CS.dRdS_i(i+1,j,:,:), CS.stable_cell(i+1,j,:),             &
%             CS.uPoL(I,j,:), CS.uPoR(I,j,:), CS.uKoL(I,j,:), CS.uKoR(I,j,:), CS.uhEff(I,j,:), &
%             CS.ppoly_coeffs_T(i,j,:,:), CS.ppoly_coeffs_S(i,j,:,:),                          &
%             CS.ppoly_coeffs_T(i+1,j,:,:), CS.ppoly_coeffs_S(i+1,j,:,:))
%       end
%     end
%   enddo ; enddo
% 
%   % Neutral surface factors at V points
%   do J = GRD.jsc-1, GRD.jec ; do i = GRD.isc, GRD.iec
%     if (GRD.mask2dCv(i,J) > 0.) then
%       if (CS.continuous_reconstruction) then
%         call find_neutral_surface_positions_continuous(GRD.ke,                                  &
%                 CS.Pint(i,j,:), CS.Tint(i,j,:), CS.Sint(i,j,:), CS.dRdT(i,j,:), CS.dRdS(i,j,:),           &
%                 CS.Pint(i,j+1,:), CS.Tint(i,j+1,:), CS.Sint(i,j+1,:), CS.dRdT(i,j+1,:), CS.dRdS(i,j+1,:), &
%                 CS.vPoL(i,J,:), CS.vPoR(i,J,:), CS.vKoL(i,J,:), CS.vKoR(i,J,:), CS.vhEff(i,J,:) )
%       else
%         call find_neutral_surface_positions_discontinuous(CS, GRD.ke, CS.ns(i,j)+CS.ns(i,j+1), &
%             CS.Pint(i,j,:), h(i,j,:), CS.T_i(i,j,:,:), CS.S_i(i,j,:,:),                      &
%             CS.dRdT_i(i,j,:,:), CS.dRdS_i(i,j,:,:), CS.stable_cell(i,j,:),                   &
%             CS.Pint(i,j+1,:), h(i,j+1,:), CS.T_i(i,j+1,:,:), CS.S_i(i,j+1,:,:),              &
%             CS.dRdT_i(i,j+1,:,:), CS.dRdS_i(i,j+1,:,:), CS.stable_cell(i,j+1,:),             &
%             CS.vPoL(I,j,:), CS.vPoR(I,j,:), CS.vKoL(I,j,:), CS.vKoR(I,j,:), CS.vhEff(I,j,:), &
%             CS.ppoly_coeffs_T(i,j,:,:), CS.ppoly_coeffs_S(i,j,:,:),                          &
%             CS.ppoly_coeffs_T(i,j+1,:,:), CS.ppoly_coeffs_S(i,j+1,:,:))
% 
%       end
%     end
%   enddo ; enddo
% 
%   % Continuous reconstructions calculate hEff as the difference between the pressures of the
%   % neutral surfaces which need to be reconverted to thickness units. The discontinuous version
%   % calculates hEff from the fraction of the nondimensional fraction of the layer occupied by
%   % the... (Please finish this thought. -RWH)
%   if (CS.continuous_reconstruction) then
%     do k = 1, CS.nsurf-1 ; do j = GRD.jsc, GRD.jec ; do I = GRD.isc-1, GRD.iec
%       if (GRD.mask2dCu(I,j) > 0.) CS.uhEff(I,j,k) = CS.uhEff(I,j,k) * pa_to_H
%     enddo ; enddo ; enddo
%     do k = 1, CS.nsurf-1 ; do J = GRD.jsc-1, GRD.jec ; do i = GRD.isc, GRD.iec
%       if (GRD.mask2dCv(i,J) > 0.) CS.vhEff(i,J,k) = CS.vhEff(i,J,k) * pa_to_H
%     enddo ; enddo ; enddo
%   end
% 
%   if (CS.id_uhEff_2d>0) then
%     hEff_sum(:,:) = 0.
%     do k = 1,CS.nsurf-1 ; do j=GRD.jsc,GRD.jec ; do i=GRD.isc-1,GRD.iec
%       hEff_sum(i,j) = hEff_sum(i,j) + CS.uhEff(i,j,k)
%     enddo ; enddo ; enddo
%     call post_data(CS.id_uhEff_2d, hEff_sum, CS.diag)
%   end
%   if (CS.id_vhEff_2d>0) then
%     hEff_sum(:,:) = 0.
%     do k = 1,CS.nsurf-1 ; do j=GRD.jsc-1,GRD.jec ; do i=GRD.isc,GRD.iec
%       hEff_sum(i,j) = hEff_sum(i,j) + CS.vhEff(i,j,k)
%     enddo ; enddo ; enddo
%     call post_data(CS.id_vhEff_2d, hEff_sum, CS.diag)
%   end
% 
% end subroutine neutral_diffusion_calc_coeffs ==========================================


% ====================================================================================
%> Update tracer concentration due to neutral diffusion; layer thickness unchanged by this update.
%% subroutine neutral_diffusion(G, GV, h, Coef_x, Coef_y, dt, Reg, CS)
  type(ocean_grid_type),                     intent(in)    :: G      %< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     %< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h      %< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: Coef_x %< dt * Kh * dy / dx at u-points [m2]
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: Coef_y %< dt * Kh * dx / dy at v-points [m2]
  real,                                      intent(in)    :: dt     %< Tracer time step * I_numitts
                                                                     %% (I_numitts in tracer_hordiff)
  type(tracer_registry_type),                pointer       :: Reg    %< Tracer registry
  type(neutral_diffusion_CS),                pointer       :: CS     %< Neutral diffusion control structure

  % Local variables
  real, dimension(SZIB_(G),SZJ_(G),CS.nsurf-1) :: uFlx        % Zonal flux of tracer [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJB_(G),CS.nsurf-1) :: vFlx        % Meridional flux of tracer
                                                              % [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJ_(G),GRD.ke)        :: tendency    % tendency array for diagn
  real, dimension(SZI_(G),SZJ_(G))             :: tendency_2d % depth integrated content tendency for diagn
  real, dimension(SZIB_(G),SZJ_(G))            :: trans_x_2d  % depth integrated diffusive tracer x-transport diagn
  real, dimension(SZI_(G),SZJB_(G))            :: trans_y_2d  % depth integrated diffusive tracer y-transport diagn
  real, dimension(GRD.ke)                        :: dTracer     % change in tracer concentration due to ndiffusion

  type(tracer_type), pointer                   :: Tracer => NULL() % Pointer to the current tracer

  integer :: i, j, k, m, ks, nk
  real :: Idt
  real :: h_neglect, h_neglect_edge

  %### Try replacing both of these with GV.H_subroundoff
  h_neglect_edge = GV.m_to_H*1.0e-10
  h_neglect = GV.m_to_H*1.0e-30

  nk = GV.ke

  do m = 1,Reg%ntr % Loop over tracer registry

    tracer => Reg%Tr(m)

    % for diagnostics
    if (tracer%id_dfxy_conc    > 0 .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 .or. &
       tracer%id_dfx_2d       > 0 .or. tracer%id_dfy_2d > 0) then
       Idt              = 1.0/dt
       tendency(:,:,:)  = 0.0
    end

    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.

    % x-flux
    do j = GRD.jsc,GRD.jec ; do I = GRD.isc-1,GRD.iec
      if (GRD.mask2dCu(I,j)>0.) then
        call neutral_surface_flux(nk, CS.nsurf, CS.deg, h(i,j,:), h(i+1,j,:),       &
                                  tracer%t(i,j,:), tracer%t(i+1,j,:), &
                                  CS.uPoL(I,j,:), CS.uPoR(I,j,:), &
                                  CS.uKoL(I,j,:), CS.uKoR(I,j,:), &
                                  CS.uhEff(I,j,:), uFlx(I,j,:), &
                                  CS.continuous_reconstruction, h_neglect, CS.remap_CS, h_neglect_edge)
      end
    enddo ; enddo

    % y-flux
    do J = GRD.jsc-1,GRD.jec ; do i = GRD.isc,GRD.iec
      if (GRD.mask2dCv(i,J)>0.) then
        call neutral_surface_flux(nk, CS.nsurf, CS.deg, h(i,j,:), h(i,j+1,:),       &
                                  tracer%t(i,j,:), tracer%t(i,j+1,:), &
                                  CS.vPoL(i,J,:), CS.vPoR(i,J,:), &
                                  CS.vKoL(i,J,:), CS.vKoR(i,J,:), &
                                  CS.vhEff(i,J,:), vFlx(i,J,:),   &
                                  CS.continuous_reconstruction, h_neglect, CS.remap_CS, h_neglect_edge)
      end
    enddo ; enddo

    % Update the tracer concentration from divergence of neutral diffusive flux components
    do j = GRD.jsc,GRD.jec ; do i = GRD.isc,GRD.iec
      if (GRD.mask2dT(i,j)>0.) then

        dTracer(:) = 0.
        do ks = 1,CS.nsurf-1
          k = CS.uKoL(I,j,ks)
          dTracer(k) = dTracer(k) + Coef_x(I,j)   * uFlx(I,j,ks)
          k = CS.uKoR(I-1,j,ks)
          dTracer(k) = dTracer(k) - Coef_x(I-1,j) * uFlx(I-1,j,ks)
          k = CS.vKoL(i,J,ks)
          dTracer(k) = dTracer(k) + Coef_y(i,J)   * vFlx(i,J,ks)
          k = CS.vKoR(i,J-1,ks)
          dTracer(k) = dTracer(k) - Coef_y(i,J-1) * vFlx(i,J-1,ks)
        enddo
        do k = 1, GV.ke
          tracer%t(i,j,k) = tracer%t(i,j,k) + dTracer(k) * &
                          ( GRD.IareaT(i,j) / ( h(i,j,k) + GV.H_subroundoff ) )
        enddo

        if (tracer%id_dfxy_conc > 0  .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 ) then
          do k = 1, GV.ke
            tendency(i,j,k) = dTracer(k) * GRD.IareaT(i,j) * Idt
          enddo
        end

      end
    enddo ; enddo

    % Diagnose vertically summed zonal flux, giving zonal tracer transport from ndiff.
    % Note sign corresponds to downgradient flux convention.
    if (tracer%id_dfx_2d > 0) then
      do j = GRD.jsc,GRD.jec ; do I = GRD.isc-1,GRD.iec
        trans_x_2d(I,j) = 0.
        if (GRD.mask2dCu(I,j)>0.) then
          do ks = 1,CS.nsurf-1
            trans_x_2d(I,j) = trans_x_2d(I,j) - Coef_x(I,j) * uFlx(I,j,ks)
          enddo
          trans_x_2d(I,j) = trans_x_2d(I,j) * Idt
        end
      enddo ; enddo
      call post_data(tracer%id_dfx_2d, trans_x_2d(:,:), CS.diag)
    end

    % Diagnose vertically summed merid flux, giving meridional tracer transport from ndiff.
    % Note sign corresponds to downgradient flux convention.
    if (tracer%id_dfy_2d > 0) then
      do J = GRD.jsc-1,GRD.jec ; do i = GRD.isc,GRD.iec
        trans_y_2d(i,J) = 0.
        if (GRD.mask2dCv(i,J)>0.) then
          do ks = 1,CS.nsurf-1
            trans_y_2d(i,J) = trans_y_2d(i,J) - Coef_y(i,J) * vFlx(i,J,ks)
          enddo
          trans_y_2d(i,J) = trans_y_2d(i,J) * Idt
        end
      enddo ; enddo
      call post_data(tracer%id_dfy_2d, trans_y_2d(:,:), CS.diag)
    end

    % post tendency of tracer content
    if (tracer%id_dfxy_cont > 0) then
      call post_data(tracer%id_dfxy_cont, tendency(:,:,:), CS.diag)
    end

    % post depth summed tendency for tracer content
    if (tracer%id_dfxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j = GRD.jsc,GRD.jec ; do i = GRD.isc,GRD.iec
        do k = 1, GV.ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_dfxy_cont_2d, tendency_2d(:,:), CS.diag)
    end

    % post tendency of tracer concentration; this step must be
    % done after posting tracer content tendency, since we alter
    % the tendency array.
    if (tracer%id_dfxy_conc > 0) then
      do k = 1, GV.ke ; do j = GRD.jsc,GRD.jec ; do i = GRD.isc,GRD.iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV.H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_dfxy_conc, tendency, CS.diag)
    end
  enddo % Loop over tracer registry

%end subroutine neutral_diffusion ==========================================



%> \namespace mom_tracer_hor_diff
%%
%% \section section_intro Introduction to the module
%%
%%    This module contains subroutines that handle horizontal
%%  diffusion (i.e., isoneutral or along layer) of tracers.
%%
%%    Each of the tracers are subject to Fickian along-coordinate
%%  diffusion if Khtr is defined and positive.  The tracer diffusion
%%  can use a suitable number of iterations to guarantee stability
%%  with an arbitrarily large time step.

%end module MOM_tracer_hor_diff


