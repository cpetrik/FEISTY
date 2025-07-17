function Tracer = sub_diff_alt(GRD,Tracer,K,ni,nj,tstep)
    % K = diffusivity in m2/s
	% tstep = time step in hours
	% nt = time steps in a day
	% dt = # seconds in nt
    
    % time step
    dt = 60.0*60.0*tstep;
	nt = (60.0*60.0*24.0) / dt;
	
    % grid size
	isd = 1;
	jsd = 2; %ignore j=1 b/c land (Antarctica)
	ied = ni;
	jed = nj;
    dxtn = GRD.dxtn;
    dyte = GRD.dyte;
    %area = GRD.dat; %all grid cells have area
    area = GRD.area; %only ocean cells, non-ocean grid cells have area=1
    mask = GRD.mask;

    dfe = zeros(ni,nj);
    dfn = zeros(ni,nj);
    gradTi = zeros(ni,nj);
    gradTj = zeros(ni,nj);
    dupwind = zeros(ni,nj);

    %% Diffusion loop
    for n = 1:nt
        % Calculate biomass gradient
        %Gradient i
        for j=jsd:jed
            for i=isd:ied
                if (i == ied)
                    gradTi(i,j) = (Tracer(isd,j) - Tracer(i,j)) *mask(i,j)*mask(isd,j);
                else
                    gradTi(i,j) = (Tracer(i+1,j) - Tracer(i,j)) *mask(i,j)*mask(i+1,j);
                end
            end
        end
        %Gradient j
        for j=jsd:jed
            for i=isd:ied
                if (j < jed)
                    gradTj(i,j) = (Tracer(i,j+1) - Tracer(i,j)) *mask(i,j)*mask(i,j+1);
                else
                    gradTj(i,j) = (Tracer(ni-i+1,j) - Tracer(i,j)) *mask(i,j)*mask(ni-i+1,j);
                end
            end
        end
        
        % Westward flux
        for j = jsd:jed
            for i = isd:ied
                % define only for ocean cells
                if (mask(i,j) > 0)
                    if (i == ied)
                        dfe(i,j) = (gradTi(isd,j)-gradTi(i,j)) .*mask(i,j).* mask(isd,j);
                    else
                        dfe(i,j) = (gradTi(i+1,j)-gradTi(i,j)) .*mask(i,j).* mask(i+1,j);
                    end
                end
            end
        end

        % Northward flux
        for j = jsd:jed
            for i = isd:ied
                % define only for ocean cells
                if (mask(i,j) > 0)
                    if (j < jed)
                        dfn(i,j) = (gradTj(i,j+1)-gradTj(i,j)) .*mask(i,j).*mask(i,j+1);
                    else
                        dfn(i,j) = (gradTj(ni-i+1,j)-gradTj(i,j)) .*mask(i,j).*mask(ni-i+1,j);
                    end
                end
            end
        end

        % Combine fluxes & Update tracers
        %This probably assumes an even grid spacing, so is wrong here
        for j = jsd:jed
            for i = isd:ied
                Tracer(i,j) = Tracer(i,j) + K .* (((dt/dyte(i,j).^2).*dfe(i,j)) + ((dt/dxtn(i,j).^2).*dfn(i,j)));
            end
        end

    end

end