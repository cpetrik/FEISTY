function Tracer = sub_NAtl_advec_vel(GRD,Tracer,uvel,vvel,iids,jids,tstep)
    % U & V = velocities in m/s
    % tstep = time step in hours
    % nt = time steps in a day
    % dt = # seconds in nt

    % time step
    dt = 60.0*60.0*tstep;
    nt = (60.0*60.0*24.0) / dt;

    % grid size
    ni = length(iids)+2;
    nj = length(jids)+2;
    gis = iids(1)-1;
    gjs = jids(1)-1;
    gie = iids(end)+1;
    gje = jids(end)+1;
    
    dxtn = GRD.dxtn(gis:gie,gjs:gje);
    dyte = GRD.dyte(gis:gie,gjs:gje);
    area = GRD.area(gis:gie,gjs:gje);
    mask = GRD.mask(gis:gie,gjs:gje);

    fe = zeros(ni,nj);
    fn = zeros(ni,nj);
    upwind = zeros(ni,nj);
    
    % start and end cells
    isd = 2;
    jsd = 2;
    ied = ni-1;
    jed = nj-1;
    
    %% Sub-daily time step loop
    for n = 1:nt
        
        % Westward flux
        for j = jsd:jed
            for i = isd:ied
                velocity = 0.5*uvel(i,j);
                upos = velocity + abs(velocity);
                uneg = velocity - abs(velocity);
                % define only for ocean cells
                if (mask(i,j) > 0)
                    fe(i,j) = dyte(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(i+1,j)) .* mask(i,j) .*mask(i+1,j);
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
                    fn(i,j) = dxtn(i,j).*(upos.*Tracer(i,j) + uneg.*Tracer(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                end
            end
        end

        % Combine fluxes
        for j = jsd:jed
            for i = isd:ied
                upwind(i,j) = mask(i,j).*(fe(i-1,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
            end
        end

        % Update tracers
        for j = jsd:jed
            for i = isd:ied
                Tracer(i,j) = Tracer(i,j) + (dt.*upwind(i,j))./area(i,j);
            end
        end

    end %time substep

end