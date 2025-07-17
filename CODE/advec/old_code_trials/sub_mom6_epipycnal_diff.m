function sub_mom6_epipycnal_diff(h, dt, Tr, K, G, num_itts)
%   G          %< ocean grid structure
%   h          %< layer thickness [H ~> m or kg m-2]
%   dt         %< time step
%   Tr         %< tracer concentration
%   K          %< horizontal diffusivity
%   num_itts   %< number of iterations (usually=1)

is = G.isc ;
ie = G.iec ;
js = G.jsc ;
je = G.jec ;
Idt = 1.0/dt;

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
Trdf2d_x = zeros(ni,nj);
Trdf2d_y = zeros(ni,nj);

for itt=1:max_itt
    
    %% Determine the fluxes through the zonal (i) faces.
    for j=js:je
        %for i=is-1:ie
        for i=is:ie %can't do (is-1), could put ni here, but then would calc it twice
            %if (G.mask2dCu(i,j) > 0.5)
            if (G.mask(i,j) > 0.5)
                % Find the acceptable range of tracer concentration around this face.
%                 Tr_min_face = min(Tr(i,j), Tr(i+1,j));
%                 Tr_max_face = max(Tr(i,j), Tr(i+1,j));
                % Include all points in diffusive pairings at this face.
                Tr_Lb = Tr(i,j);
                if (i==ie)
                    Tr_Rb = Tr(is,j);
                else
                    Tr_Rb = Tr(i+1,j);
                end
%                 Tr_La = Tr_Lb;
%                 Tr_Ra = Tr_Rb;
                
                Tr_av_L = Tr_Lb;
                Tr_av_R = Tr_Rb;
                
                h_L = h; %left thickness
                h_R = h; %right thickness
                Tr_flux = I_maxitt * K * (Tr_av_L - Tr_av_R) * ...
                    ((2.0 * h_L * h_R) / (h_L + h_R));
                
                tr_flux_conv(i,j) = tr_flux_conv(i,j) - Tr_flux;
                if (i==ie)
                    tr_flux_conv(is,j) = tr_flux_conv(is,j) + Tr_flux;
                else
                    tr_flux_conv(i+1,j) = tr_flux_conv(i+1,j) + Tr_flux;
                end
                
                Trdf2d_x(i,j) = Trdf2d_x(i,j) + Tr_flux * Idt;
            end
        end
    end % i- & j- loops over zonal faces.
    
    for j=js-1:je %can do (js-1) b/c js=2
        for i=is:ie
            %if (G.mask2dCv(i,j) > 0.5)
            if (G.mask(i,j) > 0.5)
                % Determine the fluxes through the meridional (j) faces.
                % Find the acceptable range of tracer concentration around this face.
%                 Tr_min_face = min(Tr(i,j), Tr(i,j+1));
%                 Tr_max_face = max(Tr(i,j), Tr(i,j+1));
                
                % Include all points in diffusive pairings at this face.
                Tr_Lb = Tr(i,j);
                if (j==je)
                    Tr_Rb = Tr(ni-i+1,j);
                else
                    Tr_Rb = Tr(i,j+1);
                end
%                 Tr_La = Tr_Lb;
%                 Tr_Ra = Tr_Rb;
                
                Tr_av_L = Tr_Lb;
                Tr_av_R = Tr_Rb;
                
                h_L = h;
                h_R = h;
                Tr_flux = I_maxitt * ((2.0 * h_L * h_R) / (h_L + h_R)) * ...
                    K * (Tr_av_L - Tr_av_R);
                Tr_flux_3d = Tr_flux;
                
                Trdf2d_y(i,j) = Trdf2d_y(i,j) + Tr_flux * Idt;
            end
        end
    end % i- & j- loops over meridional faces.
    
    
    for i=is:ie
        for j=js-1:je
            %if (G.mask2dCv(i,j) > 0.5)
            if (G.mask(i,j) > 0.5)
                tr_flux_conv(i,j) = tr_flux_conv(i,j) - Tr_flux_3d;
                if (j==je)
                    tr_flux_conv(ni-i+1,j) = tr_flux_conv(ni-i+1,j) + Tr_flux_3d;
                else
                    tr_flux_conv(i,j+1) = tr_flux_conv(i,j+1) + Tr_flux_3d;
                end
            end
        end
    end
    
    
    for j=js:je
        for i=is:ie
            %if ((G.mask2dT(i,j) > 0.5)
            if ((G.mask(i,j) > 0.5)
                Tr(i,j) = Tr(i,j) + tr_flux_conv(i,j) / ...
                    (h*G.area(i,j));
                tr_flux_conv(i,j) = 0.0;
            end
        end
    end
    
    
    
end % Loop over iterations

end %function