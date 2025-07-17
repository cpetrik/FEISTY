%> Update tracer concentration due to neutral diffusion; layer thickness unchanged by this update.
subroutine neutral_diffusion(G, GV, h, Coef_x, Coef_y, dt, Reg, CS)
%   type(ocean_grid_type),                     intent(in)    :: G      %< Ocean grid structure
%   type(verticalGrid_type),                   intent(in)    :: GV     %< ocean vertical grid structure
%   real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h      %< Layer thickness [H ~> m or kg m-2]
%   real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: Coef_x %< dt * Kh * dy / dx at u-points [m2]
%   real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: Coef_y %< dt * Kh * dx / dy at v-points [m2]
%   real,                                      intent(in)    :: dt     %< Tracer time step * I_numitts
%                                                                      %% (I_numitts in tracer_hordiff)
%   type(tracer_registry_type),                pointer       :: Reg    %< Tracer registry
%   type(neutral_diffusion_CS),                pointer       :: CS     %< Neutral diffusion control structure
% 
%   % Local variables
%   real, dimension(SZIB_(G),SZJ_(G),CS%nsurf-1) :: uFlx        % Zonal flux of tracer [H conc ~> m conc or conc kg m-2]
%   real, dimension(SZI_(G),SZJB_(G),CS%nsurf-1) :: vFlx        % Meridional flux of tracer
%                                                               % [H conc ~> m conc or conc kg m-2]
%   real, dimension(SZI_(G),SZJ_(G),G%ke)        :: tendency    % tendency array for diagn
%   real, dimension(SZI_(G),SZJ_(G))             :: tendency_2d % depth integrated content tendency for diagn
%   real, dimension(SZIB_(G),SZJ_(G))            :: trans_x_2d  % depth integrated diffusive tracer x-transport diagn
%   real, dimension(SZI_(G),SZJB_(G))            :: trans_y_2d  % depth integrated diffusive tracer y-transport diagn
%   real, dimension(G%ke)                        :: dTracer     % change in tracer concentration due to ndiffusion
% 
%   type(tracer_type), pointer                   :: Tracer => NULL() % Pointer to the current tracer
% 
%   integer :: i, j, k, m, ks, nk
%   real :: Idt
%   real :: h_neglect, h_neglect_edge

  %### Try replacing both of these with GV%H_subroundoff
  h_neglect_edge = GV%m_to_H*1.0e-10
  h_neglect = GV%m_to_H*1.0e-30

  nk = GV%ke

  do m = 1,Reg%ntr % Loop over tracer registry

    tracer => Reg%Tr(m)

    

    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.

    

    % Update the tracer concentration from divergence of neutral diffusive flux components
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then

        dTracer(:) = 0.
        do ks = 1,CS%nsurf-1
          k = CS%uKoL(I,j,ks)
          dTracer(k) = dTracer(k) + Coef_x(I,j)   * uFlx(I,j,ks)
          k = CS%uKoR(I-1,j,ks)
          dTracer(k) = dTracer(k) - Coef_x(I-1,j) * uFlx(I-1,j,ks)
          k = CS%vKoL(i,J,ks)
          dTracer(k) = dTracer(k) + Coef_y(i,J)   * vFlx(i,J,ks)
          k = CS%vKoR(i,J-1,ks)
          dTracer(k) = dTracer(k) - Coef_y(i,J-1) * vFlx(i,J-1,ks)
        enddo
        do k = 1, GV%ke
          tracer%t(i,j,k) = tracer%t(i,j,k) + dTracer(k) * &
                          ( G%IareaT(i,j) / ( h(i,j,k) + GV%H_subroundoff ) )
        enddo

        if (tracer%id_dfxy_conc > 0  .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 ) then
          do k = 1, GV%ke
            tendency(i,j,k) = dTracer(k) * G%IareaT(i,j) * Idt
          enddo
        endif

      endif
    enddo ; enddo

    
    
    % post depth summed tendency for tracer content
    if (tracer%id_dfxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        do k = 1, GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_dfxy_cont_2d, tendency_2d(:,:), CS%diag)
    endif

    % post tendency of tracer concentration; this step must be
    % done after posting tracer content tendency, since we alter
    % the tendency array.
    if (tracer%id_dfxy_conc > 0) then
      do k = 1, GV%ke ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_dfxy_conc, tendency, CS%diag)
    endif
  enddo % Loop over tracer registry

end subroutine neutral_diffusion