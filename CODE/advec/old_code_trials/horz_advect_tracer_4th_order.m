function horz_advect_tracer_4th_order(Adv_vel, Tracer_field)

%   type(ocean_adv_vel_type),     intent(in) :: Adv_vel
%   real, dimension(isd:,jsd:,:), intent(in) :: Tracer_field 
% 
%   real,dimension(isc:iec,jsc:jec,nk) :: horz_advect_tracer_4th_order
% 
%   integer :: i, j, k
%   integer :: im1, ip2, jm1, jp2

tmask_fourth = NaN(isc-2:iec+2,jsc-2:jec+2,nk)
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmask_fourth(i,j,k) = Grd%tmask(i,j,k)
           tmask_sixth(i,j,k)  = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo


  tracer_fourth = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_fourth(i,j,k) = Tracer_field(i,j,k)
        enddo
     enddo
  enddo
  
  flux_x = 0.0
  flux_y = 0.0

  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec
           im1   = tmask_fourth(i-1,j,k)*(i-1) + (1.0-tmask_fourth(i-1,j,k))*i
           ip2   = tmask_fourth(i+2,j,k)*(i+2) + (1.0-tmask_fourth(i+2,j,k))*(i+1)    
           flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uhrho_et(i,j,k)*  &
                (a4*(tracer_fourth(i+1,j,k)+tracer_fourth(i,j,k)) + &
                 b4*(tracer_fourth(ip2,j,k)+tracer_fourth(im1,j,k)))
        enddo
     enddo

     do j=jsc-1,jec
        do i=isc,iec
           jm1   = tmask_fourth(i,j-1,k)*(j-1) + (1.0-tmask_fourth(i,j-1,k))*j
           jp2   = tmask_fourth(i,j+2,k)*(j+2) + (1.0-tmask_fourth(i,j+2,k))*(j+1)    
           flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vhrho_nt(i,j,k)*&
                (a4*(tracer_fourth(i,j+1,k)+tracer_fourth(i,j,k)) + &
                 b4*(tracer_fourth(i,jp2,k)+tracer_fourth(i,jm1,k)))
        enddo
     enddo

  enddo