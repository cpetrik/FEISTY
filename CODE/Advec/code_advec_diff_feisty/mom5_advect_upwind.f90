! PIECES OF THE MOM5 ROUTINES

! FROM OCEAN_TRACER_ADVECT.F90
! Tracer_field = tracer concentration
! Grd%dyte = cell length on east side
! Grd%dxtn = cell length on north side
! Grd%dat = cell area
! Grd%datr = 1/cell area
! Grd%tmask = land mask
! Adv_vel%uhrho_et = eastward velocity * h * rho
! Adv_vel%vhrho_nt = northward velocity * h * rho
! Adv_vel%wrho_bt = vertical velocity * rho

! Horizontal

Tracer%wrk1(isc:iec,jsc:jec,:) =  &
-horz_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))

function horz_advect_tracer_upwind(Adv_vel, Tracer_field)
  do k=1,nk

     ! i-flux
     do j=jsc,jec
        do i=isc-1,iec
           velocity = 0.5*Adv_vel%uhrho_et(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fe(i,j)  = Grd%dyte(i,j)*(upos*Tracer_field(i,j,k) + uneg*Tracer_field(i+1,j,k)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
           flux_x(i,j,k) = fe(i,j)
        enddo
     enddo

     ! j-flux
     do j=jsc-1,jec
        do i=isc,iec
           velocity = 0.5*Adv_vel%vhrho_nt(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fn(i,j)  = Grd%dxtn(i,j)*(upos*Tracer_field(i,j,k) + uneg*Tracer_field(i,j+1,k)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
           flux_y(i,j,k) = fn(i,j)
        enddo
     enddo

    ! together
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_upwind(i,j,k) = &
            Grd%tmask(i,j,k)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
        enddo
     enddo
  enddo
end function horz_advect_tracer_upwind

do k=1,nk
   do j=jsc,jec
      do i=isc,iec
         Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
      enddo
   enddo
enddo




! Vertical
Tracer%wrk1(isc:iec,jsc:jec,:) = &
-vert_advect_tracer_upwind(Adv_vel, Tracer%field(:,:,:,taum1))

function vert_advect_tracer_upwind(Adv_vel, Tracer_field)
  ft1  = 0.0
  do k=1,nk
    kp1 = min(k+1,nk)
    do j=jsc,jec
      do i=isc,iec
        velocity = 0.5*Adv_vel%wrho_bt(i,j,k)
        wpos     = velocity + abs(velocity)
        wneg     = velocity - abs(velocity)
        ft2(i,j) = (wneg*Tracer_field(i,j,k) + wpos*Tracer_field(i,j,kp1)) &
                   *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_upwind(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        ft1(i,j) = ft2(i,j)
      enddo
    enddo
  enddo
end function vert_advect_tracer_upwind


do k=1,nk
   do j=jsc,jec
      do i=isc,iec
         Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
      enddo
   enddo
enddo



! FROM OCEAN_TRACER.F90
!#######################################################################
! <SUBROUTINE NAME="update_advection_only">
!
! <DESCRIPTION>
!
! Redo tracer updates for those that use only advection--nothing else.
! This method is useful for testing advection schemes.
!
! T_prog(n)%use_only_advection==.true. ignores all boundary forcing
! and sources, so if T_prog(n)%stf or pme, rivers, sources
! are nonzero, tracer diagnostics will spuriously indicate
! non-conservation.
!
! Assume for these tests that
!  (1) vertical advection is done fully explictly in time
!  (2) pme, rivers, stf, btf, and other sources are zero
!  (3) do not use advect_sweby_all
!
! </DESCRIPTION>
!
subroutine update_advection_only(Time, Adv_vel, Dens, Thickness, T_prog, n)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_adv_vel_type),     intent(in)    :: Adv_vel
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  integer,                      intent(in)    :: n
  integer                                     :: i, j, k
  integer                                     :: taup1, taum1

  taup1 = Time%taup1
  taum1 = Time%taum1

  T_prog(n)%th_tendency = 0.0

  call horz_advect_tracer(Time, Adv_vel, Thickness, Dens, &
       T_prog(1:num_prog_tracers), T_prog(n), n, dtime, store_flux=.FALSE.)
  call vert_advect_tracer(Time, Adv_vel, Dens, Thickness, &
       T_prog(1:num_prog_tracers), T_prog(n), n, dtime)

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           T_prog(n)%field(i,j,k,taup1) =                                    &
                (Thickness%rho_dzt(i,j,k,taum1)*T_prog(n)%field(i,j,k,taum1) &
               + dtime*T_prog(n)%th_tendency(i,j,k)) * Thickness%rho_dztr(i,j,k)
        enddo
     enddo
  enddo

end subroutine update_advection_only
! </SUBROUTINE> NAME="update_advection_only">
