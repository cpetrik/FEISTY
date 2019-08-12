function Tracer = sub_mom6_advec2_diff_vel(GRD,Tracer,K,uvel,vvel,ni,nj,tstep)
% K = diffusivity in m2/s
% uvel & vvel = velocities in m/s
% Tracer = biomass in g/m
% GRD = structure of gridfile info
% ni, nj, = grid dimensions
% tstep = time step in hours


% nt = time steps in a day
% dt = # seconds in nt
% time step
dt = 60.0*60.0*tstep;
nt = (60.0*60.0*24.0) / dt;

% grid size
isd = 1;
jsd = 2;
ied = ni;
jed = nj;
isc = 1;
jsc = 2; %ignore j=1 b/c land (Antarctica)
iec = ni;
jec = nj;

dxtn = GRD.dxtn;
dyte = GRD.dyte;
%datr = 1./GRD.dat;
area = GRD.AREA;
mask = GRD.lmask;

dfe = zeros(ni,nj);
dfn = zeros(ni,nj);
gradTi = zeros(ni,nj);
gradTj = zeros(ni,nj);
dupwind = zeros(ni,nj);

fs = zeros(ni,1);
fn = zeros(ni,1);
flux_x = zeros(ni,nj);
flux_y = zeros(ni,nj);
horz_advect = zeros(ni,nj);


%% Advection & diffusion loop
for n = 1:nt
    Tracer_field = Tracer;
    
    %Advection -----------------------------------------------
    j     = jsc-1; %starting latitude - 1 = 1st lat grid cell (land)
    fs(:) = dxtn(:,j).*vvel(:,j).*0.5.*(Tracer_field(:,j+1)+Tracer_field(:,j));
    
    for j=jsc:jec
        i = iec; %isc-1; %starting longitude - 1 = end longitude
        %fw = dyte(i,j).*uvel(i,j).*0.5.*(Tracer_field(i+1,j)+Tracer_field(i,j));
        fw = dyte(i,j).*uvel(i,j).*0.5.*(Tracer_field(isc,j)+Tracer_field(i,j));
        
        for i=isc:iec
            if (i==iec)
                fe    = dyte(i,j).*uvel(i,j).*0.5.*(Tracer_field(isc,j)+Tracer_field(i,j));
            else
                fe    = dyte(i,j).*uvel(i,j).*0.5.*(Tracer_field(i+1,j)+Tracer_field(i,j));
            end
            if (j==jec)
                fn(i) = dxtn(i,j).*vvel(i,j).*0.5.*(Tracer_field(ni-i+1,j)+Tracer_field(i,j));
            else
                fn(i) = dxtn(i,j).*vvel(i,j).*0.5.*(Tracer_field(i,j+1)+Tracer_field(i,j));
            end
            flux_x(i,j) = fe;
            flux_y(i,j) = fn(i);
            horz_advect(i,j) = mask(i,j).* (fe - fw + fn(i) - fs(i)) ./area(i,j);
            fw = fe;
        end
        
        fs(:) = fn(:);
    end
    
    %Diffusion -----------------------------------------------
    % Calculate biomass gradient
    %Gradient i
    for j=jsd:jed
        for i=isd:ied
            if (i == ied)
                gradTi(i,j) = (Tracer(isd,j) - Tracer(i,j)) ./ dyte(i,j) .*mask(i,j).*mask(isd,j);
            else
                gradTi(i,j) = (Tracer(i+1,j) - Tracer(i,j)) ./ dyte(i,j) .*mask(i,j).*mask(i+1,j);
            end
        end
    end
    %Gradient j
    for j=jsd:jed
        for i=isd:ied
            if (j < jed)
                gradTj(i,j) = (Tracer(i,j+1) - Tracer(i,j)) ./ dxtn(i,j) .*mask(i,j).*mask(i,j+1);
            else
                gradTj(i,j) = (Tracer(ni-i+1,j) - Tracer(i,j)) ./ dxtn(i,j) .*mask(i,j).*mask(ni-i+1,j);
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
    Tracer = Tracer + dt.*horz_advect - ((dt.*dupwind)./area);
    
end
