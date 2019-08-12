function Tracer = sub_advect_vel_2nd_order(GRD,Tracer,uvel,vvel,ni,nj,tstep)
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
isc = 1;
jsc = 2; %ignore j=1 b/c land (Antarctica)
iec = ni;
jec = nj;
dxtn = GRD.dxtn;
dyte = GRD.dyte;
%datr = 1./GRD.dat;
area = GRD.area;
mask = GRD.mask;

fs = zeros(ni,1);
fn = zeros(ni,1);
flux_x = zeros(ni,nj);
flux_y = zeros(ni,nj);
horz_advect = zeros(ni,nj);

%% Advection loop
for n = 1:nt
    Tracer_field = Tracer;
%     fs = zeros(ni,1);
%     fn = zeros(ni,1);
%     flux_x = zeros(ni,nj);
%     flux_y = zeros(ni,nj);
%     horz_advect = zeros(ni,nj);
    
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
            horz_advect(i,j) = mask(i,j)* (fe - fw + fn(i) - fs(i)) ./area(i,j);
            fw = fe;
        end
        
        fs(:) = fn(:);
    end
    
    Tracer = Tracer + dt.*horz_advect;
    
end
