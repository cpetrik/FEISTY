function Tracer = sub_swimgrad_advec_diff_vel(GRD,Tracer,K,uvel,vvel,nu,Q,ni,nj,tstep)
    % Move fish up gradient of food with diffusion method
    % K = diffusivity in m2/s
	% U & V = velocities in m/s
    % nu = thing to maximize/swim towards
    % Q = temp-dep swim speed
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
    dfe = zeros(ni,nj);
    dfn = zeros(ni,nj);
    sfe = zeros(ni,nj);
    sfn = zeros(ni,nj);
    gradTi = zeros(ni,nj);
    gradTj = zeros(ni,nj);
    gradKi = zeros(ni,nj);
    gradKj = zeros(ni,nj);
    upwind = zeros(ni,nj);
    dupwind = zeros(ni,nj);
    supwind = zeros(ni,nj);
    
    %% Sub-daily time loop
    for n = 1:nt
        % Calculate gradient of desired quantity (nu)
        %Gradient i
        for j=jsd:jed
            for i=isd:ied
                if (i == ied)
                    gradKi(i,j) = (nu(isd,j) - nu(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                else
                    gradKi(i,j) = (nu(i+1,j) - nu(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                end
            end
        end
        %Gradient j
        for j=jsd:jed
            for i=isd:ied
                if (j < jed)
                    gradKj(i,j) = (nu(i,j+1) - nu(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                else
                    gradKj(i,j) = (nu(ni-i+1,j) - nu(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                end
            end
        end
        gradK = (gradKi + gradKj); %.* mask;
        
        
        % Calculate biomass gradient
        %Gradient i
        for j=jsd:jed
            for i=isd:ied
                if (i == ied)
                    gradTi(i,j) = (Tracer(isd,j) - Tracer(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                else
                    gradTi(i,j) = (Tracer(i+1,j) - Tracer(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                end
            end
        end
        %Gradient j
        for j=jsd:jed
            for i=isd:ied
                if (j < jed)
                    gradTj(i,j) = (Tracer(i,j+1) - Tracer(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                else
                    gradTj(i,j) = (Tracer(ni-i+1,j) - Tracer(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                end
            end
        end
        gradT = (gradTi + gradTj); %.* mask;

        diffusiv = 0.5*K;
        kpos     = diffusiv + abs(diffusiv);
        kneg     = diffusiv - abs(diffusiv);

        % Westward flux
        for j = jsd:jed
            for i = isd:ied
                velocity = 0.5 * uvel(i,j);
                upos = velocity + abs(velocity);
                uneg = velocity - abs(velocity);
                swim = 0.5 * Q(i,j);
                wpos = swim + abs(swim);
                wneg = swim - abs(swim);

                % define only for ocean cells
                if (mask(i,j) > 0)

                    if (i == ied)
                        fe(i,j) = dyte(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(isd,j)) .* mask(i,j) .*mask(isd,j);
                        dfe(i,j) = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(isd,j)) .* mask(i,j) .* mask(isd,j);
                        sfe(i,j) = dyte(i,j).*(wpos.*gradK(i,j).*Tracer(i,j) + wneg.*gradK(isd,j).*Tracer(isd,j)) ...
                            .* mask(i,j) .* mask(isd,j);
                    else
                        fe(i,j) = dyte(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        dfe(i,j) = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        sfe(i,j) = dyte(i,j).*(wpos.*gradK(i,j).*Tracer(i,j) + wneg.*gradK(i+1,j).*Tracer(i+1,j)) ...
                            .* mask(i,j) .* mask(i+1,j);
                    end

                end
            end
        end

        % Northward flux
        for j = jsd:jed
            for i = isd:ied
                velocity = 0.5 * vvel(i,j);
                upos = velocity + abs(velocity);
                uneg = velocity - abs(velocity);
                swim = 0.5 * Q(i,j);
                wpos = swim + abs(swim);
                wneg = swim - abs(swim);

                % define only for ocean cells
                if (mask(i,j) > 0)

                    if (j < jed)
                        fn(i,j) = dxtn(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn(i,j)  = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        sfn(i,j)  = dxtn(i,j).*(wpos.*gradK(i,j).*Tracer(i,j) + wneg.*gradK(i,j+1).*Tracer(i,j+1)) ...
                            .* mask(i,j) .* mask(i,j+1);
                    else
                        fn(i,j) = dxtn(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn(i,j) = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        sfn(i,j) = dxtn(i,j).*(wpos.*gradK(i,j).*Tracer(i,j) + wneg.*gradK(ni-i+1,j).*Tracer(ni-i+1,j)) ...
                            .* mask(i,j) .* mask(ni-i+1,j);
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
                        dupwind(i,j) = mask(i,j).*(dfe(i-1,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                        supwind(i,j) = mask(i,j).*(sfe(i-1,j)-sfe(i,j)+sfn(i,j-1)-sfn(i,j));
                    else
                        upwind(i,j) = mask(i,j).*(fe(ied,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                        dupwind(i,j) = mask(i,j).*(dfe(ied,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                        supwind(i,j) = mask(i,j).*(sfe(ied,j)-sfe(i,j)+sfn(i,j-1)-sfn(i,j));
                    end
                end
            end
        end

        % Update tracers
        for j = jsd:jed
            for i = isd:ied
                Tracer(i,j) = Tracer(i,j) + (dt.*upwind(i,j))./area(i,j) ...
                    - (dt.*dupwind(i,j))./area(i,j) + (dt.*supwind(i,j))./area(i,j);
            end
        end

    end

end