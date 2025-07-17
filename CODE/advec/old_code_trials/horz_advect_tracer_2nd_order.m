function horz_advect_tracer_2nd_order(Adv_vel, Tracer_field)

  type(ocean_adv_vel_type),     intent(in) :: Adv_vel
  real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field 

  real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_2nd_order
  
  real, dimension(isd:ied) :: fs, fn
  integer                  :: i, j, k
  real                     :: fe, fw

  call mpp_clock_begin(id_clock_2nd_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
  '==>Error from ocean_tracer_advect (horz_advect_tracer_2nd_order): ocean_tracer_advect_mod not initialized')
  endif 
  
  fs(:) = 0.0; fn(:) = 0.0
  do k=1,nk 
    j     = jsc-1
    fs(:) = Grd%dxtn(:,j)*Adv_vel%vhrho_nt(:,j,k)*0.5*(Tracer_field(:,j+1,k)+Tracer_field(:,j,k))
    do j=jsc,jec
      i = isc-1
      fw = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)*0.5*(Tracer_field(i+1,j,k)+Tracer_field(i,j,k))
      do i=isc,iec
        fe    = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)*0.5*(Tracer_field(i+1,j,k)+Tracer_field(i,j,k))
        fn(i) = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k)*0.5*(Tracer_field(i,j+1,k)+Tracer_field(i,j,k))
        flux_x(i,j,k) = fe
        flux_y(i,j,k) = fn(i)
        horz_advect_tracer_2nd_order(i,j,k) = Grd%tmask(i,j,k)*(fe - fw + fn(i) - fs(i))*Grd%datr(i,j)
        fw = fe
      enddo
      fs(:) = fn(:)
    enddo
  enddo