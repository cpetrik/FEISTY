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

is = G.isc ;
ie = G.iec ;
js = G.jsc ;
je = G.jec ;
%nz = GV.ke
isd = G.isd ;
ied = G.ied ;
jsd = G.jsd ;
jed = G.jed ;
IsdB = G.IsdB ;
IedB = G.IedB ;
Idt = 1.0/dt;
%nkmb = GV.nk_rho_varies

if (num_itts <= 1)
    max_itt = 1;
    I_maxitt = 1.0;
else
    max_itt = num_itts;
    I_maxitt = 1.0 / (real(max_itt));
end


%% The tracer-specific calculations start here.

% Zero out tracer tendencies.
tr_flux_conv = zeros(ni,nj);

for itt=1:max_itt
    
    for j=js:je
        for I=is-1:ie
            if (G.mask2dCu(I,j) > 0.5)
                % Determine the fluxes through the zonal (i) faces.
                % Find the acceptable range of tracer concentration around this face.
                Tr_min_face = min(Tr(i,j), Tr(i+1,j));
                Tr_max_face = max(Tr(i,j), Tr(i+1,j));
                % Include all points in diffusive pairings at this face.
                Tr_Lb = Tr(i,j);
                Tr_Rb = Tr(i+1,j);
                Tr_La = Tr_Lb;
                Tr_Ra = Tr_Rb;
                
                Tr_Lb = Tr(i,j);
                Tr_av_L = Tr_Lb;
                
                Tr_Rb = Tr(i+1,j);
                Tr_av_R = Tr_Rb;
                
                h_L = hP_Lu(j).p(I,k); %left thickness
                h_R = hP_Ru(j).p(I,k); %right thickness
                Tr_flux = I_maxitt * khdt_epi_x(i,j) * (Tr_av_L - Tr_av_R) * ...
                    ((2.0 * h_L * h_R) / (h_L + h_R));
                
                tr_flux_conv(i,j) = tr_flux_conv(i,j) - Tr_flux;
                tr_flux_conv(i+1,j) = tr_flux_conv(i+1,j) + Tr_flux;
                
                Tr.df2d_x(I,j) = Tr.df2d_x(I,j) + Tr_flux * Idt;
            end
        end
    end % i- & j- loops over zonal faces.
    
    for J=js-1:je
        for i=is:ie
            if (G.mask2dCv(i,J) > 0.5)
                % Determine the fluxes through the meridional (j) faces.
                % Find the acceptable range of tracer concentration around this face.
                Tr_min_face = min(Tr(i,j), Tr(i,j+1));
                Tr_max_face = max(Tr(i,j), Tr(i,j+1));
                
                % Include all points in diffusive pairings at this face.
                Tr_Lb = Tr(i,j);
                Tr_Rb = Tr(i,j+1);
                Tr_La = Tr_Lb;
                Tr_Ra = Tr_Rb;
                
                Tr_av_L = Tr_Lb;
                Tr_av_R = Tr_Rb;
                
                h_L = hP_Lv(J).p(i,k);
                h_R = hP_Rv(J).p(i,k);
                Tr_flux = I_maxitt * ((2.0 * h_L * h_R) / (h_L + h_R)) * ...
                    khdt_epi_y(i,J) * (Tr_av_L - Tr_av_R);
                Tr_flux_3d = Tr_flux;
                
                Tr.df2d_y(i,J) = Tr.df2d_y(i,J) + Tr_flux * Idt;
            end
        end
    end % i- & j- loops over meridional faces.
    
    
    for i=is:ie
        for J=js-1:je
            if (G.mask2dCv(i,J) > 0.5)
                tr_flux_conv(i,j) = tr_flux_conv(i,j) - Tr_flux_3d;
                tr_flux_conv(i,j+1) = tr_flux_conv(i,j+1) + Tr_flux_3d;
            end
        end
    end
    
    
    for j=js:je
        for i=is:ie
            if ((G.mask2dT(i,j) > 0.5)
                Tr(i,j) = Tr(i,j) + tr_flux_conv(i,j) / ...
                    (h(i,j)*G.areaT(i,j));
                tr_flux_conv(i,j) = 0.0;
            end
        end
    end
    
    
    
end % Loop over iterations

% end subroutine tracer_epipycnal_ML_diff