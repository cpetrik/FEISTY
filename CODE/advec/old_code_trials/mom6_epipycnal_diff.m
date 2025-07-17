%%subroutine tracer_epipycnal_ML_diff(h, dt, Tr, ntr, khdt_epi_x, khdt_epi_y, ...
%                                     G, GV, CS, tv, num_itts)
%   type(ocean_grid_type),                    intent(inout) :: G          %< ocean grid structure
%   type(verticalGrid_type),                  intent(in)    :: GV         %< ocean vertical grid structure
%   real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h          %< layer thickness [H ~> m or kg m-2]
%   real,                                     intent(in)    :: dt         %< time step
%   type(tracer_type),                        intent(inout) :: Tr(:)      %< tracer array
%   integer,                                  intent(in)    :: ntr        %< number of tracers
%   real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: khdt_epi_x %< needs a comment
%   real, dimension(SZI_(G),SZJB_(G)),        intent(in)    :: khdt_epi_y %< needs a comment
%   type(tracer_hor_diff_CS),                 intent(inout) :: CS         %< module control structure
%   type(thermo_var_ptrs),                    intent(in)    :: tv         %< thermodynamic structure
%   integer,                                  intent(in)    :: num_itts   %< number of iterations (usually=1)
% 
% 
%   real, dimension(SZI_(G), SZJ_(G)) :: &
%     Rml_max  % The maximum coordinate density within the mixed layer [kg m-3].
%   real, dimension(SZI_(G), SZJ_(G), max(1,GV%nk_rho_varies)) :: &
%     rho_coord % The coordinate density that is used to mix along [kg m-3].
% 
%   % The naming mnemnonic is a=above,b=below,L=Left,R=Right,u=u-point,v=v-point.
%   % These are 1-D arrays of pointers to 2-d arrays to minimize memory usage.
%   type(p2d), dimension(SZJ_(G)) :: &
%     deep_wt_Lu, deep_wt_Ru, &  % The relative weighting of the deeper of a pair [nondim].
%     hP_Lu, hP_Ru       % The total thickness on each side for each pair [H ~> m or kg m-2].
% 
%   type(p2d), dimension(SZJB_(G)) :: &
%     deep_wt_Lv, deep_wt_Rv, & % The relative weighting of the deeper of a pair [nondim].
%     hP_Lv, hP_Rv       % The total thickness on each side for each pair [H ~> m or kg m-2].
% 
%   type(p2di), dimension(SZJ_(G)) :: &
%     k0b_Lu, k0a_Lu, &  % The original k-indices of the layers that participate
%     k0b_Ru, k0a_Ru     % in each pair of mixing at u-faces.
%   type(p2di), dimension(SZJB_(G)) :: &
%     k0b_Lv, k0a_Lv, &  % The original k-indices of the layers that participate
%     k0b_Rv, k0a_Rv     % in each pair of mixing at v-faces.
% 
%   real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
%     tr_flux_conv  % The flux convergence of tracers [conc H m2 ~> conc m3 or conc kg]
%   real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: Tr_flux_3d, Tr_adj_vert_L, Tr_adj_vert_R
% 
%   real, dimension(SZI_(G), SZK_(G), SZJ_(G)) :: &
%     rho_srt, & % The density of each layer of the sorted columns [kg m-3].
%     h_srt      % The thickness of each layer of the sorted columns [H ~> m or kg m-2].
%   integer, dimension(SZI_(G), SZK_(G), SZJ_(G)) :: &
%     k0_srt     % The original k-index that each layer of the sorted column
%                % corresponds to.
% 
%   real, dimension(SZK_(G)) :: &
%     h_demand_L, & % The thickness in the left (_L) or right (_R) column that
%     h_demand_R, & % is demanded to match the thickness in the counterpart [H ~> m or kg m-2].
%     h_used_L, &   % The summed thickness from the left or right columns that
%     h_used_R, &   % have actually been used [H ~> m or kg m-2].
%     h_supply_frac_L, &  % The fraction of the demanded thickness that can
%     h_supply_frac_R     % actually be supplied from a layer.
%   integer, dimension(SZK_(G)) :: &
%     kbs_Lp, &   % The sorted indicies of the Left and Right columns for
%     kbs_Rp      % each pairing.
% 
%   integer, dimension(SZI_(G), SZJ_(G))  :: &
%     num_srt, &   % The number of layers that are sorted in each column.
%     k_end_srt, & % The maximum index in each column that might need to be
%                  % sorted, based on neighboring values of max_kRho
%     max_kRho     % The index of the layer whose target density is just denser
%                  % than the densest part of the mixed layer.
%   integer, dimension(SZJ_(G))           :: &
%     max_srt      % The maximum value of num_srt in a k-row.
%   integer, dimension(SZIB_(G), SZJ_(G)) :: &
%     nPu          % The number of epipycnal pairings at each u-point.
%   integer, dimension(SZI_(G), SZJB_(G)) :: &
%     nPv          % The number of epipycnal pairings at each v-point.
%   real :: h_exclude    % A thickness that layers must attain to be considered
%                        % for inclusion in mixing [H ~> m or kg m-2].
%   real :: Idt        % The inverse of the time step [s-1].
%   real :: I_maxitt   % The inverse of the maximum number of iterations.
%   real :: rho_pair, rho_a, rho_b  % Temporary densities [kg m-3].
%   real :: Tr_min_face  % The minimum and maximum tracer concentrations
%   real :: Tr_max_face  % associated with a pairing [Conc]
%   real :: Tr_La, Tr_Lb % The 4 tracer concentrations that might be
%   real :: Tr_Ra, Tr_Rb % associated with a pairing [Conc]
%   real :: Tr_av_L    % The average tracer concentrations on the left and right
%   real :: Tr_av_R    % sides of a pairing [Conc].
%   real :: Tr_flux    % The tracer flux from left to right in a pair [conc H m2 ~> conc m3 or conc kg].
%   real :: Tr_adj_vert  % A downward vertical adjustment to Tr_flux between the
%                      % two cells that make up one side of the pairing [conc H m2 ~> conc m3 or conc kg].
%   real :: h_L, h_R   % Thicknesses to the left and right [H ~> m or kg m-2].
%   real :: wt_a, wt_b % Fractional weights of layers above and below [nondim].
%   real :: vol        % A cell volume or mass [H m2 ~> m3 or kg].
%   logical, dimension(SZK_(G)) :: &
%     left_set, &  % If true, the left or right point determines the density of
%     right_set    % of the trio.  If densities are exactly equal, both are true.
%   real :: tmp
%   real :: p_ref_cv(SZI_(G))
% 
%   integer :: k_max, k_min, k_test, itmp
%   integer :: i, j, k, k2, m, is, ie, js, je, nz, nkmb
%   integer :: isd, ied, jsd, jed, IsdB, IedB, k_size
%   integer :: kL, kR, kLa, kLb, kRa, kRb, nP, itt, ns, max_itt
%   integer :: PEmax_kRho
%   integer :: isv, iev, jsv, jev % The valid range of the indices.
% 
%   is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
%   isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
%   IsdB = G%IsdB ; IedB = G%IedB
%   Idt = 1.0/dt
%   nkmb = GV%nk_rho_varies
% 
%   if (num_itts <= 1) then
%     max_itt = 1 ; I_maxitt = 1.0
%   else
%     max_itt = num_itts ; I_maxitt = 1.0 / (real(max_itt))
%   endif
% 
%   do i=is-2,ie+2 ; p_ref_cv(i) = tv%P_Ref ; enddo
% 
%   call do_group_pass(CS%pass_t, G%Domain, clock=id_clock_pass)
%   % Determine which layers the mixed- and buffer-layers map into...
%   %$OMP parallel do default(shared)
%   do k=1,nkmb ; do j=js-2,je+2
%     call calculate_density(tv%T(:,j,k),tv%S(:,j,k), p_ref_cv, &
%                          rho_coord(:,j,k), is-2, ie-is+5, tv%eqn_of_state)
%   enddo ; enddo
% 
%   do j=js-2,je+2 ; do i=is-2,ie+2
%     Rml_max(i,j) = rho_coord(i,j,1)
%     num_srt(i,j) = 0 ; max_kRho(i,j) = 0
%   enddo ; enddo
%   do k=2,nkmb ; do j=js-2,je+2 ; do i=is-2,ie+2
%     if (Rml_max(i,j) < rho_coord(i,j,k)) Rml_max(i,j) = rho_coord(i,j,k)
%   enddo ; enddo ; enddo
%   %   Use bracketing and bisection to find the k-level that the densest of the
%   % mixed and buffer layer corresponds to, such that:
%   %     GV%Rlay(max_kRho-1) < Rml_max <= GV%Rlay(max_kRho)
% %$OMP parallel do default(none) shared(is,ie,js,je,nz,nkmb,G,GV,Rml_max,max_kRho) &
% %$OMP                          private(k_min,k_max,k_test)
%   do j=js-2,je+2 ; do i=is-2,ie+2 ; if (G%mask2dT(i,j) > 0.5) then
%     if (Rml_max(i,j) > GV%Rlay(nz)) then ; max_kRho(i,j) = nz+1
%     elseif (Rml_max(i,j) <= GV%Rlay(nkmb+1)) then ; max_kRho(i,j) = nkmb+1
%     else
%       k_min = nkmb+2 ; k_max = nz
%       do
%         k_test = (k_min + k_max) / 2
%         if (Rml_max(i,j) <= GV%Rlay(k_test-1)) then ; k_max = k_test-1
%         elseif (GV%Rlay(k_test) < Rml_max(i,j)) then ; k_min = k_test+1
%         else ; max_kRho(i,j) = k_test ; exit ; endif
% 
%         if (k_min == k_max) then ; max_kRho(i,j) = k_max ; exit ; endif
%       enddo
%     endif
%   endif ; enddo ; enddo
% 
%   PEmax_kRho = 0
%   do j=js-1,je+1 ; do i=is-1,ie+1
%     k_end_srt(i,j) = max(max_kRho(i,j), max_kRho(i-1,j), max_kRho(i+1,j), &
%                          max_kRho(i,j-1), max_kRho(i,j+1))
%     if (PEmax_kRho < k_end_srt(i,j)) PEmax_kRho = k_end_srt(i,j)
%   enddo ; enddo
%   if (PEmax_kRho > nz) PEmax_kRho = nz % PEmax_kRho could have been nz+1.
% 
%   h_exclude = 10.0*(GV%Angstrom_H + GV%H_subroundoff)
% %$OMP parallel default(none) shared(is,ie,js,je,nkmb,G,GV,h,h_exclude,num_srt,k0_srt, &
% %$OMP                               rho_srt,h_srt,PEmax_kRho,k_end_srt,rho_coord,max_srt) &
% %$OMP                       private(ns,tmp,itmp)
% %$OMP do
%   do j=js-1,je+1
%     do k=1,nkmb ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.5) then
%       if (h(i,j,k) > h_exclude) then
%         num_srt(i,j) = num_srt(i,j) + 1 ; ns = num_srt(i,j)
%         k0_srt(i,ns,j) = k
%         rho_srt(i,ns,j) = rho_coord(i,j,k)
%         h_srt(i,ns,j) = h(i,j,k)
%       endif
%     endif ; enddo ; enddo
%     do k=nkmb+1,PEmax_kRho ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.5) then
%       if ((k<=k_end_srt(i,j)) .and. (h(i,j,k) > h_exclude)) then
%         num_srt(i,j) = num_srt(i,j) + 1 ; ns = num_srt(i,j)
%         k0_srt(i,ns,j) = k
%         rho_srt(i,ns,j) = GV%Rlay(k)
%         h_srt(i,ns,j) = h(i,j,k)
%       endif
%     endif ; enddo ; enddo
%   enddo
%   % Sort each column by increasing density.  This should already be close,
%   % and the size of the arrays are small, so straight insertion is used.
% %$OMP do
%    do j=js-1,je+1; do i=is-1,ie+1
%     do k=2,num_srt(i,j) ; if (rho_srt(i,k,j) < rho_srt(i,k-1,j)) then
%       % The last segment needs to be shuffled earlier in the list.
%       do k2 = k,2,-1 ; if (rho_srt(i,k2,j) >= rho_srt(i,k2-1,j)) exit
%         itmp = k0_srt(i,k2-1,j) ; k0_srt(i,k2-1,j) = k0_srt(i,k2,j) ; k0_srt(i,k2,j) = itmp
%         tmp = rho_srt(i,k2-1,j) ; rho_srt(i,k2-1,j) = rho_srt(i,k2,j) ; rho_srt(i,k2,j) = tmp
%         tmp = h_srt(i,k2-1,j) ; h_srt(i,k2-1,j) = h_srt(i,k2,j) ; h_srt(i,k2,j) = tmp
%       enddo
%     endif ; enddo
%   enddo ; enddo
% %$OMP do
%   do j=js-1,je+1
%     max_srt(j) = 0
%     do i=is-1,ie+1 ; max_srt(j) = max(max_srt(j), num_srt(i,j)) ; enddo
%   enddo
% %$OMP end parallel
% 
%   do j=js,je
%     k_size = max(2*max_srt(j),1)
%     allocate(deep_wt_Lu(j)%p(IsdB:IedB,k_size))
%     allocate(deep_wt_Ru(j)%p(IsdB:IedB,k_size))
%     allocate(hP_Lu(j)%p(IsdB:IedB,k_size))
%     allocate(hP_Ru(j)%p(IsdB:IedB,k_size))
%     allocate(k0a_Lu(j)%p(IsdB:IedB,k_size))
%     allocate(k0a_Ru(j)%p(IsdB:IedB,k_size))
%     allocate(k0b_Lu(j)%p(IsdB:IedB,k_size))
%     allocate(k0b_Ru(j)%p(IsdB:IedB,k_size))
%   enddo
% 
% %$OMP parallel do default(none) shared(is,ie,js,je,G,num_srt,rho_srt,k0b_Lu,k0_srt, &
% %$OMP                                  k0b_Ru,k0a_Lu,k0a_Ru,deep_wt_Lu,deep_wt_Ru,  &
% %$OMP                                  h_srt,nkmb,nPu,hP_Lu,hP_Ru)                  &
% %$OMP                          private(h_demand_L,h_used_L,h_demand_R,h_used_R,     &
% %$OMP                                  kR,kL,nP,rho_pair,kbs_Lp,kbs_Rp,rho_a,rho_b, &
% %$OMP                                  wt_b,left_set,right_set,h_supply_frac_R,     &
% %$OMP                                  h_supply_frac_L)
%   do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.5) then
%     % Set up the pairings for fluxes through the zonal faces.
% 
%     do k=1,num_srt(i,j)   ; h_demand_L(k) = 0.0 ; h_used_L(k) = 0.0 ; enddo
%     do k=1,num_srt(i+1,j) ; h_demand_R(k) = 0.0 ; h_used_R(k) = 0.0 ; enddo
% 
%     % First merge the left and right lists into a single, sorted list.
% 
%     %   Discard any layers that are lighter than the lightest in the other
%     % column.  They can only participate in mixing as the lighter part of a
%     % pair of points.
%     if (rho_srt(i,1,j) < rho_srt(i+1,1,j)) then
%       kR = 1
%       do kL=2,num_srt(i,j) ; if (rho_srt(i,kL,j) >= rho_srt(i+1,1,j)) exit ; enddo
%     elseif (rho_srt(i+1,1,j) < rho_srt(i,1,j)) then
%       kL = 1
%       do kR=2,num_srt(i+1,j) ; if (rho_srt(i+1,kR,j) >= rho_srt(i,1,j)) exit ; enddo
%     else
%       kL = 1 ; kR = 1
%     endif
%     nP = 0
%     do % Loop to accumulate pairs of columns.
%       if ((kL > num_srt(i,j)) .or. (kR > num_srt(i+1,j))) exit
% 
%       if (rho_srt(i,kL,j) > rho_srt(i+1,kR,j)) then
%       % The right point is lighter and defines the density for this trio.
%         nP = nP+1 ; k = nP
%         rho_pair = rho_srt(i+1,kR,j)
% 
%         k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
%         k0a_Lu(j)%p(I,k) = k0_srt(i,kL-1,j) ; k0a_Ru(j)%p(I,k) = k0b_Ru(j)%p(I,k)
%         kbs_Lp(k) = kL ; kbs_Rp(k) = kR
% 
%         rho_a = rho_srt(i,kL-1,j) ; rho_b = rho_srt(i,kL,j)
%         wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
%           wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
%         deep_wt_Lu(j)%p(I,k) = wt_b ; deep_wt_Ru(j)%p(I,k) = 1.0
% 
%         h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j) * wt_b
%         h_demand_L(kL-1) = h_demand_L(kL-1) + 0.5*h_srt(i+1,kR,j) * (1.0-wt_b)
% 
%         kR = kR+1 ; left_set(k) = .false. ; right_set(k) = .true.
%       elseif (rho_srt(i,kL,j) < rho_srt(i+1,kR,j)) then
%       % The left point is lighter and defines the density for this trio.
%         nP = nP+1 ; k = nP
%         rho_pair = rho_srt(i,kL,j)
%         k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
%         k0a_Lu(j)%p(I,k) = k0b_Lu(j)%p(I,k) ; k0a_Ru(j)%p(I,k) = k0_srt(i+1,kR-1,j)
% 
%         kbs_Lp(k) = kL ; kbs_Rp(k) = kR
% 
%         rho_a = rho_srt(i+1,kR-1,j) ; rho_b = rho_srt(i+1,kR,j)
%         wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
%           wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
%         deep_wt_Lu(j)%p(I,k) = 1.0 ; deep_wt_Ru(j)%p(I,k) = wt_b
% 
%         h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j) * wt_b
%         h_demand_R(kR-1) = h_demand_R(kR-1) + 0.5*h_srt(i,kL,j) * (1.0-wt_b)
% 
%         kL = kL+1 ; left_set(k) = .true. ; right_set(k) = .false.
%       elseif ((k0_srt(i,kL,j) <= nkmb) .or. (k0_srt(i+1,kR,j) <= nkmb)) then
%         % The densities are exactly equal and one layer is above the interior.
%         nP = nP+1 ; k = nP
%         k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
%         k0a_Lu(j)%p(I,k) = k0b_Lu(j)%p(I,k) ; k0a_Ru(j)%p(I,k) = k0b_Ru(j)%p(I,k)
%         kbs_Lp(k) = kL ; kbs_Rp(k) = kR
%         deep_wt_Lu(j)%p(I,k) = 1.0 ; deep_wt_Ru(j)%p(I,k) = 1.0
% 
%         h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j)
%         h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)
% 
%         kL = kL+1 ; kR = kR+1 ; left_set(k) = .true. ; right_set(k) = .true.
%       else % The densities are exactly equal and in the interior.
%         % Mixing in this case has already occurred, so accumulate the thickness
%         % demanded for that mixing and skip onward.
%         h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j)
%         h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)
% 
%         kL = kL+1 ; kR = kR+1
%       endif
%     enddo % Loop to accumulate pairs of columns.
%     nPu(I,j) = nP % This is the number of active pairings.
% 
%     % Determine what fraction of the thickness "demand" can be supplied.
%     do k=1,num_srt(i+1,j)
%       h_supply_frac_R(k) = 1.0
%       if (h_demand_R(k) > 0.5*h_srt(i+1,k,j)) &
%         h_supply_frac_R(k) = 0.5*h_srt(i+1,k,j) / h_demand_R(k)
%     enddo
%     do k=1,num_srt(i,j)
%       h_supply_frac_L(k) = 1.0
%       if (h_demand_L(k) > 0.5*h_srt(i,k,j)) &
%         h_supply_frac_L(k) = 0.5*h_srt(i,k,j) / h_demand_L(k)
%     enddo
% 
%     %  Distribute the "exported" thicknesses proportionately.
%     do k=1,nPu(I,j)
%       kL = kbs_Lp(k) ; kR = kbs_Rp(k)
%       hP_Lu(j)%p(I,k) = 0.0 ; hP_Ru(j)%p(I,k) = 0.0
%       if (left_set(k)) then % Add the contributing thicknesses on the right.
%         if (deep_wt_Ru(j)%p(I,k) < 1.0) then
%           hP_Ru(j)%p(I,k) = 0.5*h_srt(i,kL,j) * min(h_supply_frac_R(kR), h_supply_frac_R(kR-1))
%           wt_b = deep_wt_Ru(j)%p(I,k)
%           h_used_R(kR-1) = h_used_R(kR-1) + (1.0 - wt_b)*hP_Ru(j)%p(I,k)
%           h_used_R(kR) = h_used_R(kR) + wt_b*hP_Ru(j)%p(I,k)
%         else
%           hP_Ru(j)%p(I,k) = 0.5*h_srt(i,kL,j) * h_supply_frac_R(kR)
%           h_used_R(kR) = h_used_R(kR) + hP_Ru(j)%p(I,k)
%         endif
%       endif
%       if (right_set(k)) then % Add the contributing thicknesses on the left.
%         if (deep_wt_Lu(j)%p(I,k) < 1.0) then
%           hP_Lu(j)%p(I,k) = 0.5*h_srt(i+1,kR,j) * min(h_supply_frac_L(kL), h_supply_frac_L(kL-1))
%           wt_b = deep_wt_Lu(j)%p(I,k)
%           h_used_L(kL-1) = h_used_L(kL-1) + (1.0 - wt_b)*hP_Lu(j)%p(I,k)
%           h_used_L(kL) = h_used_L(kL) + wt_b*hP_Lu(j)%p(I,k)
%         else
%           hP_Lu(j)%p(I,k) = 0.5*h_srt(i+1,kR,j) * h_supply_frac_L(kL)
%           h_used_L(kL) = h_used_L(kL) + hP_Lu(j)%p(I,k)
%         endif
%       endif
%     enddo
% 
%     %   The left-over thickness (at least half the layer thickness) is now
%     % added to the thicknesses of the importing columns.
%     do k=1,nPu(I,j)
%       if (left_set(k)) hP_Lu(j)%p(I,k) = hP_Lu(j)%p(I,k) + &
%                            (h_srt(i,kbs_Lp(k),j) - h_used_L(kbs_Lp(k)))
%       if (right_set(k)) hP_Ru(j)%p(I,k) = hP_Ru(j)%p(I,k) + &
%                             (h_srt(i+1,kbs_Rp(k),j) - h_used_R(kbs_Rp(k)))
%     enddo
% 
%   endif ; enddo ; enddo % i- & j- loops over zonal faces.
% 
%   do J=js-1,je
%     k_size = max(max_srt(j)+max_srt(j+1),1)
%     allocate(deep_wt_Lv(J)%p(isd:ied,k_size))
%     allocate(deep_wt_Rv(J)%p(isd:ied,k_size))
%     allocate(hP_Lv(J)%p(isd:ied,k_size))
%     allocate(hP_Rv(J)%p(isd:ied,k_size))
%     allocate(k0a_Lv(J)%p(isd:ied,k_size))
%     allocate(k0a_Rv(J)%p(isd:ied,k_size))
%     allocate(k0b_Lv(J)%p(isd:ied,k_size))
%     allocate(k0b_Rv(J)%p(isd:ied,k_size))
%   enddo
% 
% %$OMP parallel do default(none) shared(is,ie,js,je,G,num_srt,rho_srt,k0b_Lv,k0b_Rv, &
% %$OMP                                  k0_srt,k0a_Lv,k0a_Rv,deep_wt_Lv,deep_wt_Rv,  &
% %$OMP                                  h_srt,nkmb,nPv,hP_Lv,hP_Rv)                  &
% %$OMP                          private(h_demand_L,h_used_L,h_demand_R,h_used_R,     &
% %$OMP                                  kR,kL,nP,rho_pair,kbs_Lp,kbs_Rp,rho_a,rho_b, &
% %$OMP                                  wt_b,left_set,right_set,h_supply_frac_R,     &
% %$OMP                                  h_supply_frac_L)
%   do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.5) then
%     % Set up the pairings for fluxes through the meridional faces.
% 
%     do k=1,num_srt(i,j)   ; h_demand_L(k) = 0.0 ; h_used_L(k) = 0.0 ; enddo
%     do k=1,num_srt(i,j+1) ; h_demand_R(k) = 0.0 ; h_used_R(k) = 0.0 ; enddo
% 
%     % First merge the left and right lists into a single, sorted list.
% 
%     %   Discard any layers that are lighter than the lightest in the other
%     % column.  They can only participate in mixing as the lighter part of a
%     % pair of points.
%     if (rho_srt(i,1,j) < rho_srt(i,1,j+1)) then
%       kR = 1
%       do kL=2,num_srt(i,j) ; if (rho_srt(i,kL,j) >= rho_srt(i,1,j+1)) exit ; enddo
%     elseif (rho_srt(i,1,j+1) < rho_srt(i,1,j)) then
%       kL = 1
%       do kR=2,num_srt(i,j+1) ; if (rho_srt(i,kR,j+1) >= rho_srt(i,1,j)) exit ; enddo
%     else
%       kL = 1 ; kR = 1
%     endif
%     nP = 0
%     do % Loop to accumulate pairs of columns.
%       if ((kL > num_srt(i,j)) .or. (kR > num_srt(i,j+1))) exit
% 
%       if (rho_srt(i,kL,j) > rho_srt(i,kR,j+1)) then
%       % The right point is lighter and defines the density for this trio.
%         nP = nP+1 ; k = nP
%         rho_pair = rho_srt(i,kR,j+1)
% 
%         k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j)   ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
%         k0a_Lv(J)%p(i,k) = k0_srt(i,kL-1,j) ; k0a_Rv(J)%p(i,k) = k0b_Rv(J)%p(i,k)
%         kbs_Lp(k) = kL ; kbs_Rp(k) = kR
% 
%         rho_a = rho_srt(i,kL-1,j) ; rho_b = rho_srt(i,kL,j)
%         wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
%           wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
%         deep_wt_Lv(J)%p(i,k) = wt_b ; deep_wt_Rv(J)%p(i,k) = 1.0
% 
%         h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1) * wt_b
%         h_demand_L(kL-1) = h_demand_L(kL-1) + 0.5*h_srt(i,kR,j+1) * (1.0-wt_b)
% 
%         kR = kR+1 ; left_set(k) = .false. ; right_set(k) = .true.
%       elseif (rho_srt(i,kL,j) < rho_srt(i,kR,j+1)) then
%       % The left point is lighter and defines the density for this trio.
%         nP = nP+1 ; k = nP
%         rho_pair = rho_srt(i,kL,j)
%         k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j) ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
%         k0a_Lv(J)%p(i,k) = k0b_Lv(J)%p(i,k) ; k0a_Rv(J)%p(i,k) = k0_srt(i,kR-1,j+1)
% 
%         kbs_Lp(k) = kL ; kbs_Rp(k) = kR
% 
%         rho_a = rho_srt(i,kR-1,j+1) ; rho_b = rho_srt(i,kR,j+1)
%         wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
%           wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
%         deep_wt_Lv(J)%p(i,k) = 1.0 ; deep_wt_Rv(J)%p(i,k) = wt_b
% 
%         h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j) * wt_b
%         h_demand_R(kR-1) = h_demand_R(kR-1) + 0.5*h_srt(i,kL,j) * (1.0-wt_b)
% 
%         kL = kL+1 ; left_set(k) = .true. ; right_set(k) = .false.
%       elseif ((k0_srt(i,kL,j) <= nkmb) .or. (k0_srt(i,kR,j+1) <= nkmb)) then
%         % The densities are exactly equal and one layer is above the interior.
%         nP = nP+1 ; k = nP
%         k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j) ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
%         k0a_Lv(J)%p(i,k) = k0b_Lv(J)%p(i,k)  ; k0a_Rv(J)%p(i,k) = k0b_Rv(J)%p(i,k)
%         kbs_Lp(k) = kL ; kbs_Rp(k) = kR
%         deep_wt_Lv(J)%p(i,k) = 1.0 ; deep_wt_Rv(J)%p(i,k) = 1.0
% 
%         h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1)
%         h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)
% 
%         kL = kL+1 ; kR = kR+1 ; left_set(k) = .true. ; right_set(k) = .true.
%       else % The densities are exactly equal and in the interior.
%         % Mixing in this case has already occurred, so accumulate the thickness
%         % demanded for that mixing and skip onward.
%         h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1)
%         h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)
% 
%         kL = kL+1 ; kR = kR+1
%       endif
%     enddo % Loop to accumulate pairs of columns.
%     nPv(i,J) = nP % This is the number of active pairings.
% 
%     % Determine what fraction of the thickness "demand" can be supplied.
%     do k=1,num_srt(i,j+1)
%       h_supply_frac_R(k) = 1.0
%       if (h_demand_R(k) > 0.5*h_srt(i,k,j+1)) &
%         h_supply_frac_R(k) = 0.5*h_srt(i,k,j+1) / h_demand_R(k)
%     enddo
%     do k=1,num_srt(i,j)
%       h_supply_frac_L(k) = 1.0
%       if (h_demand_L(k) > 0.5*h_srt(i,k,j)) &
%         h_supply_frac_L(k) = 0.5*h_srt(i,k,j) / h_demand_L(k)
%     enddo
% 
%     %  Distribute the "exported" thicknesses proportionately.
%     do k=1,nPv(i,J)
%       kL = kbs_Lp(k) ; kR = kbs_Rp(k)
%       hP_Lv(J)%p(i,k) = 0.0 ; hP_Rv(J)%p(i,k) = 0.0
%       if (left_set(k)) then % Add the contributing thicknesses on the right.
%         if (deep_wt_Rv(J)%p(i,k) < 1.0) then
%           hP_Rv(J)%p(i,k) = 0.5*h_srt(i,kL,j) * min(h_supply_frac_R(kR), h_supply_frac_R(kR-1))
%           wt_b = deep_wt_Rv(J)%p(i,k)
%           h_used_R(kR-1) = h_used_R(kR-1) + (1.0 - wt_b) * hP_Rv(J)%p(i,k)
%           h_used_R(kR) = h_used_R(kR) + wt_b * hP_Rv(J)%p(i,k)
%         else
%           hP_Rv(J)%p(i,k) = 0.5*h_srt(i,kL,j) * h_supply_frac_R(kR)
%           h_used_R(kR) = h_used_R(kR) + hP_Rv(J)%p(i,k)
%         endif
%       endif
%       if (right_set(k)) then % Add the contributing thicknesses on the left.
%         if (deep_wt_Lv(J)%p(i,k) < 1.0) then
%           hP_Lv(J)%p(i,k) = 0.5*h_srt(i,kR,j+1) * min(h_supply_frac_L(kL), h_supply_frac_L(kL-1))
%           wt_b = deep_wt_Lv(J)%p(i,k)
%           h_used_L(kL-1) = h_used_L(kL-1) + (1.0 - wt_b) * hP_Lv(J)%p(i,k)
%           h_used_L(kL) = h_used_L(kL) + wt_b * hP_Lv(J)%p(i,k)
%         else
%           hP_Lv(J)%p(i,k) = 0.5*h_srt(i,kR,j+1) * h_supply_frac_L(kL)
%           h_used_L(kL) = h_used_L(kL) + hP_Lv(J)%p(i,k)
%         endif
%       endif
%     enddo
% 
%     %   The left-over thickness (at least half the layer thickness) is now
%     % added to the thicknesses of the importing columns.
%     do k=1,nPv(i,J)
%       if (left_set(k)) hP_Lv(J)%p(i,k) = hP_Lv(J)%p(i,k) + &
%                             (h_srt(i,kbs_Lp(k),j) - h_used_L(kbs_Lp(k)))
%       if (right_set(k)) hP_Rv(J)%p(i,k) = hP_Rv(J)%p(i,k) + &
%                              (h_srt(i,kbs_Rp(k),j+1) - h_used_R(kbs_Rp(k)))
%     enddo
% 
% 
%   endif ; enddo ; enddo % i- & j- loops over meridional faces.

%% The tracer-specific calculations start here.

  % Zero out tracer tendencies.
  do k=1,PEmax_kRho ; do j=js-1,je+1 ; do i=is-1,ie+1
    tr_flux_conv(i,j,k) = 0.0
  enddo ; enddo ; enddo

  do itt=1,max_itt

    if (itt > 1) then % The halos have already been filled if itt==1.
      call do_group_pass(CS%pass_t, G%Domain, clock=id_clock_pass)
    endif

    do m=1,ntr
%$OMP parallel do default(none) shared(is,ie,js,je,G,Tr,nkmb,nPu,m,max_kRho,nz,h,h_exclude, &
%$OMP                                  k0b_Lu,k0b_Ru,deep_wt_Lu,k0a_Lu,deep_wt_Ru,k0a_Ru,   &
%$OMP                                  hP_Lu,hP_Ru,I_maxitt,khdt_epi_x,tr_flux_conv,Idt) &
%$OMP                          private(Tr_min_face,Tr_max_face,kLa,kLb,kRa,kRb,Tr_La, &
%$OMP                                     Tr_Lb,Tr_Ra,Tr_Rb,Tr_av_L,wt_b,Tr_av_R,h_L,h_R, &
%$OMP                                     Tr_flux,Tr_adj_vert,wt_a,vol)
      do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.5) then
        % Determine the fluxes through the zonal faces.

        % Find the acceptable range of tracer concentration around this face.
        if (nPu(I,j) >= 1) then
          Tr_min_face = min(Tr(m).t(i,j,1), Tr(m).t(i+1,j,1))
          Tr_max_face = max(Tr(m).t(i,j,1), Tr(m).t(i+1,j,1))
          do k=2,nkmb
            Tr_min_face = min(Tr_min_face, Tr(m).t(i,j,k), Tr(m).t(i+1,j,k))
            Tr_max_face = max(Tr_max_face, Tr(m).t(i,j,k), Tr(m).t(i+1,j,k))
          enddo

          % Include the next two layers denser than the densest buffer layer.
          kLa = nkmb+1 ; if (max_kRho(i,j) < nz+1) kLa = max_kRho(i,j)
          kLb = kLa ; if (max_kRho(i,j) < nz) kLb = max_kRho(i,j)+1
          kRa = nkmb+1 ; if (max_kRho(i+1,j) < nz+1) kRa = max_kRho(i+1,j)
          kRb = kRa ; if (max_kRho(i+1,j) < nz) kRb = max_kRho(i+1,j)+1
          Tr_La = Tr_min_face ; Tr_Lb = Tr_La ; Tr_Ra = Tr_La ; Tr_Rb = Tr_La
          if (h(i,j,kLa) > h_exclude) Tr_La = Tr(m).t(i,j,kLa)
          if (h(i,j,kLb) > h_exclude) Tr_La = Tr(m).t(i,j,kLb)
          if (h(i+1,j,kRa) > h_exclude) Tr_Ra = Tr(m).t(i+1,j,kRa)
          if (h(i+1,j,kRb) > h_exclude) Tr_Rb = Tr(m).t(i+1,j,kRb)
          Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)

          % Include all points in diffusive pairings at this face.
          do k=1,nPu(I,j)
            Tr_Lb = Tr(m).t(i,j,k0b_Lu(j)%p(I,k))
            Tr_Rb = Tr(m).t(i+1,j,k0b_Ru(j)%p(I,k))
            Tr_La = Tr_Lb ; Tr_Ra = Tr_Rb
            if (deep_wt_Lu(j)%p(I,k) < 1.0) Tr_La = Tr(m).t(i,j,k0a_Lu(j)%p(I,k))
            if (deep_wt_Ru(j)%p(I,k) < 1.0) Tr_Ra = Tr(m).t(i+1,j,k0a_Ru(j)%p(I,k))
            Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
            Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          enddo
        endif

        do k=1,nPu(I,j)
          kLb = k0b_Lu(j)%p(I,k) ; Tr_Lb = Tr(m).t(i,j,kLb) ; Tr_av_L = Tr_Lb
          if (deep_wt_Lu(j)%p(I,k) < 1.0) then
            kLa = k0a_Lu(j)%p(I,k) ; Tr_La = Tr(m).t(i,j,kLa)
            wt_b = deep_wt_Lu(j)%p(I,k)
            Tr_av_L = wt_b*Tr_Lb + (1.0-wt_b)*Tr_La
          endif

          kRb = k0b_Ru(j)%p(I,k) ; Tr_Rb = Tr(m).t(i+1,j,kRb) ; Tr_av_R = Tr_Rb
          if (deep_wt_Ru(j)%p(I,k) < 1.0) then
            kRa = k0a_Ru(j)%p(I,k) ; Tr_Ra = Tr(m).t(i+1,j,kRa)
            wt_b = deep_wt_Ru(j)%p(I,k)
            Tr_av_R = wt_b*Tr_Rb + (1.0-wt_b)*Tr_Ra
          endif

          h_L = hP_Lu(j)%p(I,k) ; h_R = hP_Ru(j)%p(I,k)
          Tr_flux = I_maxitt * khdt_epi_x(I,j) * (Tr_av_L - Tr_av_R) * &
            ((2.0 * h_L * h_R) / (h_L + h_R))


          if (deep_wt_Lu(j)%p(I,k) >= 1.0) then
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Lu(j)%p(I,k) ; wt_a = 1.0 - wt_b
            vol = hP_Lu(j)%p(I,k) * G%areaT(i,j)

            %   Ensure that the tracer flux does not drive the tracer values
            % outside of the range Tr_min_face <= Tr <= Tr_max_face, or if it
            % does that the concentration in both contributing peices exceed
            % this range equally. With downgradient fluxes and the initial tracer
            % concentrations determining the valid range, the latter condition
            % only enters for large values of the effective diffusive CFL number.
            if (Tr_flux > 0.0) then
              if (Tr_La < Tr_Lb) then ; if (vol*(Tr_La-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = -wt_a * min(Tr_flux - vol * (Tr_La-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Lb - Tr_La))
              else ; if (vol*(Tr_Lb-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = wt_b * min(Tr_flux - vol * (Tr_Lb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_La - Tr_Lb))
              endif
            elseif (Tr_flux < 0.0) then
              if (Tr_La > Tr_Lb) then ; if (vol * (Tr_max_face-Tr_La) < -Tr_flux) &
                Tr_adj_vert = wt_a * min(-Tr_flux - vol * (Tr_max_face-Tr_La), &
                                         (vol*wt_b) * (Tr_La - Tr_Lb))
              else ; if (vol*(Tr_max_face-Tr_Lb) < -Tr_flux) &
                Tr_adj_vert = -wt_b * min(-Tr_flux - vol * (Tr_max_face-Tr_Lb), &
                                          (vol*wt_a)*(Tr_Lb - Tr_La))
              endif
            endif

            tr_flux_conv(i,j,kLa) = tr_flux_conv(i,j,kLa) - (wt_a*Tr_flux + Tr_adj_vert)
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - (wt_b*Tr_flux - Tr_adj_vert)
          endif

          if (deep_wt_Ru(j)%p(I,k) >= 1.0) then
            tr_flux_conv(i+1,j,kRb) = tr_flux_conv(i+1,j,kRb) + Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Ru(j)%p(I,k) ; wt_a = 1.0 - wt_b
            vol = hP_Ru(j)%p(I,k) * G%areaT(i+1,j)

            %   Ensure that the tracer flux does not drive the tracer values
            % outside of the range Tr_min_face <= Tr <= Tr_max_face, or if it
            % does that the concentration in both contributing peices exceed
            % this range equally. With downgradient fluxes and the initial tracer
            % concentrations determining the valid range, the latter condition
            % only enters for large values of the effective diffusive CFL number.
            if (Tr_flux < 0.0) then
              if (Tr_Ra < Tr_Rb) then ; if (vol * (Tr_Ra-Tr_min_face) < -Tr_flux) &
                Tr_adj_vert = -wt_a * min(-Tr_flux - vol * (Tr_Ra-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Rb - Tr_Ra))
              else ; if (vol*(Tr_Rb-Tr_min_face) < (-Tr_flux)) &
                Tr_adj_vert = wt_b * min(-Tr_flux - vol * (Tr_Rb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_Ra - Tr_Rb))
              endif
            elseif (Tr_flux > 0.0) then
              if (Tr_Ra > Tr_Rb) then ; if (vol * (Tr_max_face-Tr_Ra) < Tr_flux) &
                Tr_adj_vert = wt_a * min(Tr_flux - vol * (Tr_max_face-Tr_Ra), &
                                         (vol*wt_b) * (Tr_Ra - Tr_Rb))
              else ; if (vol*(Tr_max_face-Tr_Rb) < Tr_flux) &
                Tr_adj_vert = -wt_b * min(Tr_flux - vol * (Tr_max_face-Tr_Rb), &
                                          (vol*wt_a)*(Tr_Rb - Tr_Ra))
              endif
            endif

            tr_flux_conv(i+1,j,kRa) = tr_flux_conv(i+1,j,kRa) + &
                                            (wt_a*Tr_flux - Tr_adj_vert)
            tr_flux_conv(i+1,j,kRb) = tr_flux_conv(i+1,j,kRb) + &
                                            (wt_b*Tr_flux + Tr_adj_vert)
          endif
          if (associated(Tr(m).df2d_x)) &
            Tr(m).df2d_x(I,j) = Tr(m).df2d_x(I,j) + Tr_flux * Idt
        enddo % Loop over pairings at faces.
      endif ; enddo ; enddo % i- & j- loops over zonal faces.

%$OMP parallel do default(none) shared(is,ie,js,je,G,Tr,nkmb,nPv,m,max_kRho,nz,h,h_exclude, &
%$OMP                                  k0b_Lv,k0b_Rv,deep_wt_Lv,k0a_Lv,deep_wt_Rv,k0a_Rv,   &
%$OMP                                  hP_Lv,hP_Rv,I_maxitt,khdt_epi_y,Tr_flux_3d,          &
%$OMP                                  Tr_adj_vert_L,Tr_adj_vert_R,Idt)                     &
%$OMP                          private(Tr_min_face,Tr_max_face,kLa,kLb,kRa,kRb,             &
%$OMP                                  Tr_La,Tr_Lb,Tr_Ra,Tr_Rb,Tr_av_L,wt_b,Tr_av_R,        &
%$OMP                                  h_L,h_R,Tr_flux,Tr_adj_vert,wt_a,vol)
      do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.5) then
        % Determine the fluxes through the meridional faces.

        % Find the acceptable range of tracer concentration around this face.
        if (nPv(i,J) >= 1) then
          Tr_min_face = min(Tr(m).t(i,j,1), Tr(m).t(i,j+1,1))
          Tr_max_face = max(Tr(m).t(i,j,1), Tr(m).t(i,j+1,1))
          do k=2,nkmb
            Tr_min_face = min(Tr_min_face, Tr(m).t(i,j,k), Tr(m).t(i,j+1,k))
            Tr_max_face = max(Tr_max_face, Tr(m).t(i,j,k), Tr(m).t(i,j+1,k))
          enddo

          % Include the next two layers denser than the densest buffer layer.
          kLa = nkmb+1 ; if (max_kRho(i,j) < nz+1) kLa = max_kRho(i,j)
          kLb = kLa ; if (max_kRho(i,j) < nz) kLb = max_kRho(i,j)+1
          kRa = nkmb+1 ; if (max_kRho(i,j+1) < nz+1) kRa = max_kRho(i,j+1)
          kRb = kRa ; if (max_kRho(i,j+1) < nz) kRb = max_kRho(i,j+1)+1
          Tr_La = Tr_min_face ; Tr_Lb = Tr_La ; Tr_Ra = Tr_La ; Tr_Rb = Tr_La
          if (h(i,j,kLa) > h_exclude) Tr_La = Tr(m).t(i,j,kLa)
          if (h(i,j,kLb) > h_exclude) Tr_La = Tr(m).t(i,j,kLb)
          if (h(i,j+1,kRa) > h_exclude) Tr_Ra = Tr(m).t(i,j+1,kRa)
          if (h(i,j+1,kRb) > h_exclude) Tr_Rb = Tr(m).t(i,j+1,kRb)
          Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)

          % Include all points in diffusive pairings at this face.
          do k=1,nPv(i,J)
            Tr_Lb = Tr(m).t(i,j,k0b_Lv(J)%p(i,k)) ; Tr_Rb = Tr(m).t(i,j+1,k0b_Rv(J)%p(i,k))
            Tr_La = Tr_Lb ; Tr_Ra = Tr_Rb
            if (deep_wt_Lv(J)%p(i,k) < 1.0) Tr_La = Tr(m).t(i,j,k0a_Lv(J)%p(i,k))
            if (deep_wt_Rv(J)%p(i,k) < 1.0) Tr_Ra = Tr(m).t(i,j+1,k0a_Rv(J)%p(i,k))
            Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
            Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          enddo
        endif

        do k=1,nPv(i,J)
          kLb = k0b_Lv(J)%p(i,k) ; Tr_Lb = Tr(m).t(i,j,kLb) ; Tr_av_L = Tr_Lb
          if (deep_wt_Lv(J)%p(i,k) < 1.0) then
            kLa = k0a_Lv(J)%p(i,k) ; Tr_La = Tr(m).t(i,j,kLa)
            wt_b = deep_wt_Lv(J)%p(i,k)
            Tr_av_L = wt_b * Tr_Lb + (1.0-wt_b) * Tr_La
          endif

          kRb = k0b_Rv(J)%p(i,k) ; Tr_Rb = Tr(m).t(i,j+1,kRb) ; Tr_av_R = Tr_Rb
          if (deep_wt_Rv(J)%p(i,k) < 1.0) then
            kRa = k0a_Rv(J)%p(i,k) ; Tr_Ra = Tr(m).t(i,j+1,kRa)
            wt_b = deep_wt_Rv(J)%p(i,k)
            Tr_av_R = wt_b * Tr_Rb + (1.0-wt_b) * Tr_Ra
          endif

          h_L = hP_Lv(J)%p(i,k) ; h_R = hP_Rv(J)%p(i,k)
          Tr_flux = I_maxitt * ((2.0 * h_L * h_R) / (h_L + h_R)) * &
                    khdt_epi_y(i,J) * (Tr_av_L - Tr_av_R)
          Tr_flux_3d(i,j,k) = Tr_flux

          if (deep_wt_Lv(J)%p(i,k) < 1.0) then
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Lv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            vol = hP_Lv(J)%p(i,k) * G%areaT(i,j)

            %   Ensure that the tracer flux does not drive the tracer values
            % outside of the range Tr_min_face <= Tr <= Tr_max_face.
            if (Tr_flux > 0.0) then
              if (Tr_La < Tr_Lb) then ; if (vol * (Tr_La-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = -wt_a * min(Tr_flux - vol * (Tr_La-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Lb - Tr_La))
              else ; if (vol*(Tr_Lb-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = wt_b * min(Tr_flux - vol * (Tr_Lb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_La - Tr_Lb))
              endif
            elseif (Tr_flux < 0.0) then
              if (Tr_La > Tr_Lb) then ; if (vol * (Tr_max_face-Tr_La) < -Tr_flux) &
                Tr_adj_vert = wt_a * min(-Tr_flux - vol * (Tr_max_face-Tr_La), &
                                         (vol*wt_b) * (Tr_La - Tr_Lb))
              else ; if (vol*(Tr_max_face-Tr_Lb) < -Tr_flux) &
                Tr_adj_vert = -wt_b * min(-Tr_flux - vol * (Tr_max_face-Tr_Lb), &
                                          (vol*wt_a)*(Tr_Lb - Tr_La))
              endif
            endif
            Tr_adj_vert_L(i,j,k) = Tr_adj_vert
          endif

          if (deep_wt_Rv(J)%p(i,k) < 1.0) then
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Rv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            vol = hP_Rv(J)%p(i,k) * G%areaT(i,j+1)

            %   Ensure that the tracer flux does not drive the tracer values
            % outside of the range Tr_min_face <= Tr <= Tr_max_face.
            if (Tr_flux < 0.0) then
              if (Tr_Ra < Tr_Rb) then ; if (vol * (Tr_Ra-Tr_min_face) < -Tr_flux) &
                Tr_adj_vert = -wt_a * min(-Tr_flux - vol * (Tr_Ra-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Rb - Tr_Ra))
              else ; if (vol*(Tr_Rb-Tr_min_face) < (-Tr_flux)) &
                Tr_adj_vert = wt_b * min(-Tr_flux - vol * (Tr_Rb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_Ra - Tr_Rb))
              endif
            elseif (Tr_flux > 0.0) then
              if (Tr_Ra > Tr_Rb) then ; if (vol * (Tr_max_face-Tr_Ra) < Tr_flux) &
                Tr_adj_vert = wt_a * min(Tr_flux - vol * (Tr_max_face-Tr_Ra), &
                                         (vol*wt_b) * (Tr_Ra - Tr_Rb))
              else ; if (vol*(Tr_max_face-Tr_Rb) < Tr_flux) &
                Tr_adj_vert = -wt_b * min(Tr_flux - vol * (Tr_max_face-Tr_Rb), &
                                          (vol*wt_a)*(Tr_Rb - Tr_Ra))
              endif
            endif
            Tr_adj_vert_R(i,j,k) = Tr_adj_vert
          endif
          if (associated(Tr(m).df2d_y)) &
            Tr(m).df2d_y(i,J) = Tr(m).df2d_y(i,J) + Tr_flux * Idt
        enddo % Loop over pairings at faces.
      endif ; enddo ; enddo % i- & j- loops over meridional faces.
%$OMP parallel do default(none) shared(is,ie,js,je,G,nPv,k0b_Lv,k0b_Rv,deep_wt_Lv,  &
%$OMP                                  tr_flux_conv,Tr_flux_3d,k0a_Lv,Tr_adj_vert_L,&
%$OMP                                  deep_wt_Rv,k0a_Rv,Tr_adj_vert_R) &
%$OMP                          private(kLa,kLb,kRa,kRb,wt_b,wt_a)
      do i=is,ie ; do J=js-1,je ; if (G%mask2dCv(i,J) > 0.5) then
        do k=1,nPv(i,J)
          kLb = k0b_Lv(J)%p(i,k); kRb = k0b_Rv(J)%p(i,k)
          if (deep_wt_Lv(J)%p(i,k) >= 1.0) then
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - Tr_flux_3d(i,j,k)
          else
            kLa = k0a_Lv(J)%p(i,k)
            wt_b = deep_wt_Lv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            tr_flux_conv(i,j,kLa) = tr_flux_conv(i,j,kLa) - (wt_a*Tr_flux_3d(i,j,k) + Tr_adj_vert_L(i,j,k))
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - (wt_b*Tr_flux_3d(i,j,k) - Tr_adj_vert_L(i,j,k))
          endif
          if (deep_wt_Rv(J)%p(i,k) >= 1.0) then
            tr_flux_conv(i,j+1,kRb) = tr_flux_conv(i,j+1,kRb) + tr_flux_3d(i,j,k)
          else
            kRa = k0a_Rv(J)%p(i,k)
            wt_b = deep_wt_Rv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            tr_flux_conv(i,j+1,kRa) = tr_flux_conv(i,j+1,kRa) + &
                                            (wt_a*Tr_flux_3d(i,j,k) - Tr_adj_vert_R(i,j,k))
            tr_flux_conv(i,j+1,kRb) = tr_flux_conv(i,j+1,kRb) + &
                                            (wt_b*Tr_flux_3d(i,j,k) + Tr_adj_vert_R(i,j,k))
          endif
        enddo
      endif ; enddo ; enddo
%$OMP parallel do default(none) shared(PEmax_kRho,is,ie,js,je,G,h,Tr,tr_flux_conv,m)
      do k=1,PEmax_kRho ; do j=js,je ; do i=is,ie
        if ((G%mask2dT(i,j) > 0.5) .and. (h(i,j,k) > 0.0)) then
          Tr(m).t(i,j,k) = Tr(m).t(i,j,k) + tr_flux_conv(i,j,k) / &
                                            (h(i,j,k)*G%areaT(i,j))
          tr_flux_conv(i,j,k) = 0.0
        endif
      enddo ; enddo ; enddo

    enddo % Loop over tracers
  enddo % Loop over iterations

%   do j=js,je
%     deallocate(deep_wt_Lu(j)%p) ; deallocate(deep_wt_Ru(j)%p)
%     deallocate(Hp_Lu(j)%p)  ; deallocate(Hp_Ru(j)%p)
%     deallocate(k0a_Lu(j)%p) ; deallocate(k0a_Ru(j)%p)
%     deallocate(k0b_Lu(j)%p) ; deallocate(k0b_Ru(j)%p)
%   enddo
% 
%   do J=js-1,je
%     deallocate(deep_wt_Lv(J)%p) ; deallocate(deep_wt_Rv(J)%p)
%     deallocate(Hp_Lv(J)%p)  ; deallocate(Hp_Rv(J)%p)
%     deallocate(k0a_Lv(J)%p) ; deallocate(k0a_Rv(J)%p)
%     deallocate(k0b_Lv(J)%p) ; deallocate(k0b_Rv(J)%p)
%   enddo

% end subroutine tracer_epipycnal_ML_diff