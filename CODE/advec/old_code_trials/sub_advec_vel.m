function Tracer = sub_advec_vel(GRD,Tracer,uvel,vvel,ni,nj,tstep)
    % K = diffusivity in m2/s
	% U & V = velocities in m/s
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
    %area = GRD.dat;
    area = GRD.area;
    mask = GRD.mask;

    fe = zeros(ni,nj);
    fn = zeros(ni,nj);
    upwind = zeros(ni,nj);

    %% Advection loop
    for n = 1:nt
        % Westward flux
        for j = jsd:jed
            for i = isd:ied
                velocity = 0.5*uvel(i,j);
                upos = velocity + abs(velocity);
                uneg = velocity - abs(velocity);
                % define only for ocean cells
                if (mask(i,j) > 0)
                    if (i == ied)
                        fe(i,j) = dyte(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(isd,j)) .* mask(i,j) .*mask(isd,j);
                    else
                        fe(i,j) = dyte(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                    end
                end
            end
        end

        % Northward flux
        for j = jsd:jed
            for i = isd:ied
                velocity = 0.5*vvel(i,j);
                upos = velocity + abs(velocity);
                uneg = velocity - abs(velocity);
                % define only for ocean cells
                if (mask(i,j) > 0)
                    if (j < jed)
                        fn(i,j) = dxtn(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                    else
                        fn(i,j) = dxtn(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                    end
                end
            end
        end

        % Combine fluxes
        for j = jsd:jed
            for i = isd:ied
                if (j > 1)
                    if (i > 1)
                        upwind(i,j) = mask(i,j).*(fe(i-1,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                    else
                        upwind(i,j) = mask(i,j).*(fe(ied,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                    end
                end
            end
        end

        % Update tracers
        for j = jsd:jed
            for i = isd:ied
                Tracer(i,j) = Tracer(i,j) + (dt.*upwind(i,j))./area(i,j);
            end
        end

    end

end