function Tracer = sub_advect_vel_4th_order(GRD,Tracer,K,uvel,vvel,ni,nj,tstep)

% U & V = velocities in m/s
% tstep = time step in hours
% nt = time steps in a day
% dt = # seconds in nt

% constants
a4 = 7/12;
b4 = -1/12;

% time step
dt = 60.0*60.0*tstep;
nt = (60.0*60.0*24.0) / dt;

% grid size
isd = 1;
jsd = 1;
ied = ni;
jed = nj;
%add 2 cell buffer to everything, so everything shifted up 2
isc = 3;
jsc = 3;
iec = ni+2;
jec = nj+2;
dxtn = GRD.dxtn;
dyte = GRD.dyte;
%area = GRD.dat;
area = GRD.area;
datr = 1./area;
mask = GRD.mask;


%% add 2 cell buffer/wrapper to everything
%make sure the value there reflects i+1 when i=ied is isd, etc.
advect = zeros(ni+4,nj+4);

% allocate tmask_fourth(isc-2:iec+2,jsc-2:jec+2)
tmask_fourth = zeros(iec+2,jec+2);
datr4 = zeros(iec+2,jec+2);
Uvel = zeros(iec+2,jec+2);
Vvel = zeros(iec+2,jec+2);
dx4 = zeros(iec+2,jec+2);
dy4 = zeros(iec+2,jec+2);
for j=jsc:jec
    for i=isc:iec
        tmask_fourth(i,j) = mask(i-2,j-2);
        datr4(i,j) = datr(i-2,j-2);
        Uvel(i,j) = uvel(i-2,j-2);
        Vvel(i,j) = vvel(i-2,j-2);
        dx4(i,j) = dxtn(i-2,j-2);
        dy4(i,j) = dyte(i-2,j-2);
    end
end

%START
%when i=1, really ied-1
tmask_fourth(1,jsc:jec) = mask(ied-1,:);
datr4(1,jsc:jec)      = datr(ied-1,:);
Uvel(1,jsc:jec)       = uvel(ied-1,:);
Vvel(1,jsc:jec)       = vvel(ied-1,:);
dx4(1,jsc:jec)        = dxtn(ied-1,:);
dy4(1,jsc:jec)        = dyte(ied-1,:);
%when i=2, really ied
tmask_fourth(2,jsc:jec) = mask(ied,:);
datr4(2,jsc:jec)      = datr(ied,:);
Uvel(2,jsc:jec)       = uvel(ied,:);
Vvel(2,jsc:jec)       = vvel(ied,:);
dx4(2,jsc:jec)        = dxtn(ied,:);
dy4(2,jsc:jec)        = dyte(ied,:);
%when i=3, really isd

%when j=1, really jed-1
%when j=2, really jed
%when j=3, really jsd
%The first 4 latitudes are all Antarctic land, so keep at zero

%END
%when i=iec=ni+2, really ied
%when i=ni+3, really isd
tmask_fourth(ni+3,jsc:jec) = mask(isd,:);
datr4(ni+3,jsc:jec)      = datr(isd,:);
Uvel(ni+3,jsc:jec)       = uvel(isd,:);
Vvel(ni+3,jsc:jec)       = vvel(isd,:);
dx4(ni+3,jsc:jec)        = dxtn(isd,:);
dy4(ni+3,jsc:jec)        = dyte(isd,:);
%when i=ni+4, really isd+1
tmask_fourth(ni+4,jsc:jec) = mask(isd+1,:);
datr4(ni+4,jsc:jec)      = datr(isd+1,:);
Uvel(ni+4,jsc:jec)       = uvel(isd+1,:);
Vvel(ni+4,jsc:jec)       = vvel(isd+1,:);
dx4(ni+4,jsc:jec)        = dxtn(isd+1,:);
dy4(ni+4,jsc:jec)        = dyte(isd+1,:);

%when j=jec=nj+2, really jed
%when j=nj+3, really (i,j+1) ==> (ni-i+1,j)  %cross N pole
for i=isd:ied
    tmask_fourth(i+2,nj+3) = mask(ni-i+1,jed);
    datr4(i+2,nj+3) = datr(ni-i+1,jed);
    Uvel(i+2,nj+3) = uvel(ni-i+1,jed);
    Vvel(i+2,nj+3) = vvel(ni-i+1,jed);
    dx4(i+2,nj+3) = dxtn(ni-i+1,jed);
    dy4(i+2,nj+3) = dyte(ni-i+1,jed);
end
%when j=nj+4, really (i,j+2) ==> (ni-i+1,j-1) %wrap over N pole
for i=isd:ied
    tmask_fourth(i+2,nj+4) = mask(ni-i+1,jed-1);
    datr4(i+2,nj+4) = datr(ni-i+1,jed-1);
    Uvel(i+2,nj+4) = uvel(ni-i+1,jed-1);
    Vvel(i+2,nj+4) = vvel(ni-i+1,jed-1);
    dx4(i+2,nj+4) = dxtn(ni-i+1,jed-1);
    dy4(i+2,nj+4) = dyte(ni-i+1,jed-1);
end

%% Advection loop
for n = 1:nt
    Tracer_field = Tracer;
    
    %Advection -----------------------------------------------
    tracer_fourth = zeros(iec+2,jec+2);
    for j=jsc:jec
        for i=isc:iec
            tracer_fourth(i,j) = Tracer_field(i-2,j-2);
        end
    end
    
    flux_x = zeros(iec+2,jec+2);
    flux_y = zeros(iec+2,jec+2);
    
    for j=jsc:jec
        for i=isc-1:iec
            im1 = tmask_fourth(i-1,j)*(i-1) + (1.0-tmask_fourth(i-1,j))*i;
            ip2 = tmask_fourth(i+2,j)*(i+2) + (1.0-tmask_fourth(i+2,j))*(i+1);
            flux_x(i,j) = dy4(i,j).*Uvel(i,j).*  ...
                (a4*(tracer_fourth(i+1,j)+tracer_fourth(i,j)) + ...
                b4*(tracer_fourth(ip2,j)+tracer_fourth(im1,j)));
        end
    end
    
    for j=jsc-1:jec
        for i=isc:iec
            jm1 = tmask_fourth(i,j-1)*(j-1) + (1.0-tmask_fourth(i,j-1))*j;
            jp2 = tmask_fourth(i,j+2)*(j+2) + (1.0-tmask_fourth(i,j+2))*(j+1);
            flux_y(i,j) = dx4(i,j).*Vvel(i,j).* ...
                (a4*(tracer_fourth(i,j+1)+tracer_fourth(i,j)) + ...
                b4*(tracer_fourth(i,jp2)+tracer_fourth(i,jm1)));
        end
    end
    
    for j=jsc:jec
        for i=isc:iec
            advect(i,j) = tmask_fourth(i,j)* ...
                (flux_x(i,j) - flux_x(i-1,j) + flux_y(i,j) - flux_y(i,j-1))*datr4(i,j);
        end
    end
    
    
    %Diffusion -----------------------------------------------
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
                % define only for ocean cells
                if (mask(i,j) > 0)
                    if (i == ied)
                        dfe(i,j) = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(isd,j)) .* mask(i,j) .* mask(isd,j);
                    else
                        dfe(i,j) = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i+1,j)) .* mask(i,j) .* mask(i+1,j);
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
                        dfn(i,j)  = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                    else
                        dfn(i,j) = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                    end
                end
            end
        end

        % Combine fluxes
        for j = jsd:jed
            for i = isd:ied
                if (j > 1)
                    if (i > 1)
                        dupwind(i,j) = mask(i,j).*(dfe(i-1,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                    else
                        dupwind(i,j) = mask(i,j).*(dfe(ied,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                    end
                end
            end
        end

        
    % Update tracers    
    Tracer = Tracer + dt.*advect(isc:iec,jsc:jec) - ((dt.*dupwind)./area);
end


end