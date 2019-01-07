function [B1,B2,B3,B4,B5,B6,B7,B8] = sub_advec_diff_vel(GRD,K,U,V,ni,nj,tstep,B1,B2,B3,B4,B5,B6,B7,B8)
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
    area = GRD.area;
    mask = GRD.mask;
    id = find(mask == 1);
    
    % put vectors into COBALT grid
    uvel = zeros(ni,nj);
    vvel = zeros(ni,nj);
    uvel(id) = U;
    vvel(id) = V;
    T1 = zeros(ni,nj);
    T2 = zeros(ni,nj);
    T3 = zeros(ni,nj);
    T4 = zeros(ni,nj);
    T5 = zeros(ni,nj);
    T6 = zeros(ni,nj);
    T7 = zeros(ni,nj);
    T8 = zeros(ni,nj);
    T1(id) = B1;
    T2(id) = B2;
    T3(id) = B3;
    T4(id) = B4;
    T5(id) = B5;
    T6(id) = B6;
    T7(id) = B7;
    T8(id) = B8;
    
    fe1 = zeros(ni,nj);
    fe2 = zeros(ni,nj);
    fe3 = zeros(ni,nj);
    fe4 = zeros(ni,nj);
    fe5 = zeros(ni,nj);
    fe6 = zeros(ni,nj);
    fe7 = zeros(ni,nj);
    fe8 = zeros(ni,nj);
    fn1 = zeros(ni,nj);
    fn2 = zeros(ni,nj);
    fn3 = zeros(ni,nj);
    fn4 = zeros(ni,nj);
    fn5 = zeros(ni,nj);
    fn6 = zeros(ni,nj);
    fn7 = zeros(ni,nj);
    fn8 = zeros(ni,nj);
    dfe1 = zeros(ni,nj);
    dfe2 = zeros(ni,nj);
    dfe3 = zeros(ni,nj);
    dfe4 = zeros(ni,nj);
    dfe5 = zeros(ni,nj);
    dfe6 = zeros(ni,nj);
    dfe7 = zeros(ni,nj);
    dfe8 = zeros(ni,nj);
    dfn1 = zeros(ni,nj);
    dfn2 = zeros(ni,nj);
    dfn3 = zeros(ni,nj);
    dfn4 = zeros(ni,nj);
    dfn5 = zeros(ni,nj);
    dfn6 = zeros(ni,nj);
    dfn7 = zeros(ni,nj);
    dfn8 = zeros(ni,nj);
    gradTi1 = zeros(ni,nj);
    gradTi2 = zeros(ni,nj);
    gradTi3 = zeros(ni,nj);
    gradTi4 = zeros(ni,nj);
    gradTi5 = zeros(ni,nj);
    gradTi6 = zeros(ni,nj);
    gradTi7 = zeros(ni,nj);
    gradTi8 = zeros(ni,nj);
    gradTj1 = zeros(ni,nj);
    gradTj2 = zeros(ni,nj);
    gradTj3 = zeros(ni,nj);
    gradTj4 = zeros(ni,nj);
    gradTj5 = zeros(ni,nj);
    gradTj6 = zeros(ni,nj);
    gradTj7 = zeros(ni,nj);
    gradTj8 = zeros(ni,nj);
    upwind1 = zeros(ni,nj);
    upwind2 = zeros(ni,nj);
    upwind3 = zeros(ni,nj);
    upwind4 = zeros(ni,nj);
    upwind5 = zeros(ni,nj);
    upwind6 = zeros(ni,nj);
    upwind7 = zeros(ni,nj);
    upwind8 = zeros(ni,nj);
    dupwind1 = zeros(ni,nj);
    dupwind2 = zeros(ni,nj);
    dupwind3 = zeros(ni,nj);
    dupwind4 = zeros(ni,nj);
    dupwind5 = zeros(ni,nj);
    dupwind6 = zeros(ni,nj);
    dupwind7 = zeros(ni,nj);
    dupwind8 = zeros(ni,nj);
    
    %% Advection loop
    for n = 1:nt
        % Calculate biomass gradient
        %Gradient i
        for j=jsd:jed
            for i=isd:ied
                if (i == ied)
                    gradTi1(i,j) = (T1(isd,j) - T1(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi2(i,j) = (T2(isd,j) - T2(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi3(i,j) = (T3(isd,j) - T3(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi4(i,j) = (T4(isd,j) - T4(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi5(i,j) = (T5(isd,j) - T5(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi6(i,j) = (T6(isd,j) - T6(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi7(i,j) = (T7(isd,j) - T7(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                    gradTi8(i,j) = (T8(isd,j) - T8(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
                else
                    gradTi1(i,j) = (T1(i+1,j) - T1(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi2(i,j) = (T2(i+1,j) - T2(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi3(i,j) = (T3(i+1,j) - T3(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi4(i,j) = (T4(i+1,j) - T4(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi5(i,j) = (T5(i+1,j) - T5(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi6(i,j) = (T6(i+1,j) - T6(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi7(i,j) = (T7(i+1,j) - T7(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                    gradTi8(i,j) = (T8(i+1,j) - T8(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
                end
            end
        end
        %Gradient j
        for j=jsd:jed
            for i=isd:ied
                if (j < jed)
                    gradTj1(i,j) = (T1(i,j+1) - T1(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj2(i,j) = (T2(i,j+1) - T2(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj3(i,j) = (T3(i,j+1) - T3(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj4(i,j) = (T4(i,j+1) - T4(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj5(i,j) = (T5(i,j+1) - T5(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj6(i,j) = (T6(i,j+1) - T6(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj7(i,j) = (T7(i,j+1) - T7(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                    gradTj8(i,j) = (T8(i,j+1) - T8(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
                else
                    gradTj1(i,j) = (T1(ni-i+1,j) - T1(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj2(i,j) = (T2(ni-i+1,j) - T2(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj3(i,j) = (T3(ni-i+1,j) - T3(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj4(i,j) = (T4(ni-i+1,j) - T4(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj5(i,j) = (T5(ni-i+1,j) - T5(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj6(i,j) = (T6(ni-i+1,j) - T6(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj7(i,j) = (T7(ni-i+1,j) - T7(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                    gradTj8(i,j) = (T8(ni-i+1,j) - T8(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
                end
            end
        end
        gradT1 = (gradTi1 + gradTj1); 
        gradT2 = (gradTi2 + gradTj2);
        gradT3 = (gradTi3 + gradTj3);
        gradT4 = (gradTi4 + gradTj4); 
        gradT5 = (gradTi5 + gradTj5);
        gradT6 = (gradTi6 + gradTj6);
        gradT7 = (gradTi7 + gradTj7); 
        gradT8 = (gradTi8 + gradTj8);

        diffusiv = 0.5*K;
        kpos     = diffusiv + abs(diffusiv);
        kneg     = diffusiv - abs(diffusiv);

        % Westward flux
        for j = jsd:jed
            for i = isd:ied
                velocity = 0.5*uvel(i,j);
                upos = velocity + abs(velocity);
                uneg = velocity - abs(velocity);
                % define only for ocean cells
                if (mask(i,j) > 0)
                    if (i == ied)
                        fe1(i,j) = dyte(i,j).*(upos.*T1(i,j) + uneg.*T1(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe2(i,j) = dyte(i,j).*(upos.*T2(i,j) + uneg.*T2(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe3(i,j) = dyte(i,j).*(upos.*T3(i,j) + uneg.*T3(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe4(i,j) = dyte(i,j).*(upos.*T4(i,j) + uneg.*T4(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe5(i,j) = dyte(i,j).*(upos.*T5(i,j) + uneg.*T5(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe6(i,j) = dyte(i,j).*(upos.*T6(i,j) + uneg.*T6(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe7(i,j) = dyte(i,j).*(upos.*T7(i,j) + uneg.*T7(isd,j)) .* mask(i,j) .*mask(isd,j);
                        fe8(i,j) = dyte(i,j).*(upos.*T8(i,j) + uneg.*T8(isd,j)) .* mask(i,j) .*mask(isd,j);
                        
                        dfe1(i,j) = dyte(i,j).*(kpos.*gradT1(i,j) + kneg.*gradT1(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe2(i,j) = dyte(i,j).*(kpos.*gradT2(i,j) + kneg.*gradT2(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe3(i,j) = dyte(i,j).*(kpos.*gradT3(i,j) + kneg.*gradT3(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe4(i,j) = dyte(i,j).*(kpos.*gradT4(i,j) + kneg.*gradT4(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe5(i,j) = dyte(i,j).*(kpos.*gradT5(i,j) + kneg.*gradT5(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe6(i,j) = dyte(i,j).*(kpos.*gradT6(i,j) + kneg.*gradT6(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe7(i,j) = dyte(i,j).*(kpos.*gradT7(i,j) + kneg.*gradT7(isd,j)) .* mask(i,j) .* mask(isd,j);
                        dfe8(i,j) = dyte(i,j).*(kpos.*gradT8(i,j) + kneg.*gradT8(isd,j)) .* mask(i,j) .* mask(isd,j);
                    else
                        fe1(i,j) = dyte(i,j).*(upos.*T1(i,j) + uneg.*T1(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe2(i,j) = dyte(i,j).*(upos.*T2(i,j) + uneg.*T2(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe3(i,j) = dyte(i,j).*(upos.*T3(i,j) + uneg.*T3(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe4(i,j) = dyte(i,j).*(upos.*T4(i,j) + uneg.*T4(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe5(i,j) = dyte(i,j).*(upos.*T5(i,j) + uneg.*T5(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe6(i,j) = dyte(i,j).*(upos.*T6(i,j) + uneg.*T6(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe7(i,j) = dyte(i,j).*(upos.*T7(i,j) + uneg.*T7(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        fe8(i,j) = dyte(i,j).*(upos.*T8(i,j) + uneg.*T8(i+1,j)) .* mask(i,j) .*mask(i+1,j);
                        
                        dfe1(i,j) = dyte(i,j).*(kpos.*gradT1(i,j) + kneg.*gradT1(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe2(i,j) = dyte(i,j).*(kpos.*gradT2(i,j) + kneg.*gradT2(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe3(i,j) = dyte(i,j).*(kpos.*gradT3(i,j) + kneg.*gradT3(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe4(i,j) = dyte(i,j).*(kpos.*gradT4(i,j) + kneg.*gradT4(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe5(i,j) = dyte(i,j).*(kpos.*gradT5(i,j) + kneg.*gradT5(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe6(i,j) = dyte(i,j).*(kpos.*gradT6(i,j) + kneg.*gradT6(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe7(i,j) = dyte(i,j).*(kpos.*gradT7(i,j) + kneg.*gradT7(i+1,j)) .* mask(i,j) .* mask(i+1,j);
                        dfe8(i,j) = dyte(i,j).*(kpos.*gradT8(i,j) + kneg.*gradT8(i+1,j)) .* mask(i,j) .* mask(i+1,j);
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
                        fn1(i,j) = dxtn(i,j).*(upos.*T1(i,j) + uneg.*T1(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn2(i,j) = dxtn(i,j).*(upos.*T2(i,j) + uneg.*T2(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn3(i,j) = dxtn(i,j).*(upos.*T3(i,j) + uneg.*T3(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn4(i,j) = dxtn(i,j).*(upos.*T4(i,j) + uneg.*T4(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn5(i,j) = dxtn(i,j).*(upos.*T5(i,j) + uneg.*T5(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn6(i,j) = dxtn(i,j).*(upos.*T6(i,j) + uneg.*T6(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn7(i,j) = dxtn(i,j).*(upos.*T7(i,j) + uneg.*T7(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        fn8(i,j) = dxtn(i,j).*(upos.*T8(i,j) + uneg.*T8(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        
                        dfn1(i,j)  = dxtn(i,j).*(kpos.*gradT1(i,j) + kneg.*gradT1(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn2(i,j)  = dxtn(i,j).*(kpos.*gradT2(i,j) + kneg.*gradT2(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn3(i,j)  = dxtn(i,j).*(kpos.*gradT3(i,j) + kneg.*gradT3(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn4(i,j)  = dxtn(i,j).*(kpos.*gradT4(i,j) + kneg.*gradT4(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn5(i,j)  = dxtn(i,j).*(kpos.*gradT5(i,j) + kneg.*gradT5(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn6(i,j)  = dxtn(i,j).*(kpos.*gradT6(i,j) + kneg.*gradT6(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn7(i,j)  = dxtn(i,j).*(kpos.*gradT7(i,j) + kneg.*gradT7(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                        dfn8(i,j)  = dxtn(i,j).*(kpos.*gradT8(i,j) + kneg.*gradT8(i,j+1)) .* mask(i,j) .* mask(i,j+1);
                    else
                        fn1(i,j) = dxtn(i,j).*(upos.*T1(i,j) + uneg.*T1(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn2(i,j) = dxtn(i,j).*(upos.*T2(i,j) + uneg.*T2(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn3(i,j) = dxtn(i,j).*(upos.*T3(i,j) + uneg.*T3(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn4(i,j) = dxtn(i,j).*(upos.*T4(i,j) + uneg.*T4(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn5(i,j) = dxtn(i,j).*(upos.*T5(i,j) + uneg.*T5(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn6(i,j) = dxtn(i,j).*(upos.*T6(i,j) + uneg.*T6(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn7(i,j) = dxtn(i,j).*(upos.*T7(i,j) + uneg.*T7(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        fn8(i,j) = dxtn(i,j).*(upos.*T8(i,j) + uneg.*T8(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        
                        dfn1(i,j) = dxtn(i,j).*(kpos.*gradT1(i,j) + kneg.*gradT1(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn2(i,j) = dxtn(i,j).*(kpos.*gradT2(i,j) + kneg.*gradT2(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn3(i,j) = dxtn(i,j).*(kpos.*gradT3(i,j) + kneg.*gradT3(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn4(i,j) = dxtn(i,j).*(kpos.*gradT4(i,j) + kneg.*gradT4(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn5(i,j) = dxtn(i,j).*(kpos.*gradT5(i,j) + kneg.*gradT5(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn6(i,j) = dxtn(i,j).*(kpos.*gradT6(i,j) + kneg.*gradT6(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn7(i,j) = dxtn(i,j).*(kpos.*gradT7(i,j) + kneg.*gradT7(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                        dfn8(i,j) = dxtn(i,j).*(kpos.*gradT8(i,j) + kneg.*gradT8(ni-i+1,j)) .* mask(i,j) .* mask(ni-i+1,j);
                    end
                end
            end
        end

        % Combine fluxes
        for j = jsd:jed
            for i = isd:ied
                if (j > 1)
                    if (i > 1)
                        upwind1(i,j) = mask(i,j).*(fe1(i-1,j)-fe1(i,j)+fn1(i,j-1)-fn1(i,j));
                        upwind2(i,j) = mask(i,j).*(fe2(i-1,j)-fe2(i,j)+fn2(i,j-1)-fn2(i,j));
                        upwind3(i,j) = mask(i,j).*(fe3(i-1,j)-fe3(i,j)+fn3(i,j-1)-fn3(i,j));
                        upwind4(i,j) = mask(i,j).*(fe4(i-1,j)-fe4(i,j)+fn4(i,j-1)-fn4(i,j));
                        upwind5(i,j) = mask(i,j).*(fe5(i-1,j)-fe5(i,j)+fn5(i,j-1)-fn5(i,j));
                        upwind6(i,j) = mask(i,j).*(fe6(i-1,j)-fe6(i,j)+fn6(i,j-1)-fn6(i,j));
                        upwind7(i,j) = mask(i,j).*(fe7(i-1,j)-fe7(i,j)+fn7(i,j-1)-fn7(i,j));
                        upwind8(i,j) = mask(i,j).*(fe8(i-1,j)-fe8(i,j)+fn8(i,j-1)-fn8(i,j));
                        
                        dupwind1(i,j) = mask(i,j).*(dfe1(i-1,j)-dfe1(i,j)+dfn1(i,j-1)-dfn1(i,j));
                        dupwind2(i,j) = mask(i,j).*(dfe2(i-1,j)-dfe2(i,j)+dfn2(i,j-1)-dfn2(i,j));
                        dupwind3(i,j) = mask(i,j).*(dfe3(i-1,j)-dfe3(i,j)+dfn3(i,j-1)-dfn3(i,j));
                        dupwind4(i,j) = mask(i,j).*(dfe4(i-1,j)-dfe4(i,j)+dfn4(i,j-1)-dfn4(i,j));
                        dupwind5(i,j) = mask(i,j).*(dfe5(i-1,j)-dfe5(i,j)+dfn5(i,j-1)-dfn5(i,j));
                        dupwind6(i,j) = mask(i,j).*(dfe6(i-1,j)-dfe6(i,j)+dfn6(i,j-1)-dfn6(i,j));
                        dupwind7(i,j) = mask(i,j).*(dfe7(i-1,j)-dfe7(i,j)+dfn7(i,j-1)-dfn7(i,j));
                        dupwind8(i,j) = mask(i,j).*(dfe8(i-1,j)-dfe8(i,j)+dfn8(i,j-1)-dfn8(i,j));
                    else
                        upwind1(i,j) = mask(i,j).*(fe1(ied,j)-fe1(i,j)+fn1(i,j-1)-fn1(i,j));
                        upwind2(i,j) = mask(i,j).*(fe2(ied,j)-fe2(i,j)+fn2(i,j-1)-fn2(i,j));
                        upwind3(i,j) = mask(i,j).*(fe3(ied,j)-fe3(i,j)+fn3(i,j-1)-fn3(i,j));
                        upwind4(i,j) = mask(i,j).*(fe4(ied,j)-fe4(i,j)+fn4(i,j-1)-fn4(i,j));
                        upwind5(i,j) = mask(i,j).*(fe5(ied,j)-fe5(i,j)+fn5(i,j-1)-fn5(i,j));
                        upwind6(i,j) = mask(i,j).*(fe6(ied,j)-fe6(i,j)+fn6(i,j-1)-fn6(i,j));
                        upwind7(i,j) = mask(i,j).*(fe7(ied,j)-fe7(i,j)+fn7(i,j-1)-fn7(i,j));
                        upwind8(i,j) = mask(i,j).*(fe8(ied,j)-fe8(i,j)+fn8(i,j-1)-fn8(i,j));
                        
                        dupwind1(i,j) = mask(i,j).*(dfe1(ied,j)-dfe1(i,j)+dfn1(i,j-1)-dfn1(i,j));
                        dupwind2(i,j) = mask(i,j).*(dfe2(ied,j)-dfe2(i,j)+dfn2(i,j-1)-dfn2(i,j));
                        dupwind3(i,j) = mask(i,j).*(dfe3(ied,j)-dfe3(i,j)+dfn3(i,j-1)-dfn3(i,j));
                        dupwind4(i,j) = mask(i,j).*(dfe4(ied,j)-dfe4(i,j)+dfn4(i,j-1)-dfn4(i,j));
                        dupwind5(i,j) = mask(i,j).*(dfe5(ied,j)-dfe5(i,j)+dfn5(i,j-1)-dfn5(i,j));
                        dupwind6(i,j) = mask(i,j).*(dfe6(ied,j)-dfe6(i,j)+dfn6(i,j-1)-dfn6(i,j));
                        dupwind7(i,j) = mask(i,j).*(dfe7(ied,j)-dfe7(i,j)+dfn7(i,j-1)-dfn7(i,j));
                        dupwind8(i,j) = mask(i,j).*(dfe8(ied,j)-dfe8(i,j)+dfn8(i,j-1)-dfn8(i,j));
                    end
                end
            end
        end

        % Update tracers
        for j = jsd:jed
            for i = isd:ied
                T1(i,j) = T1(i,j) + (dt.*upwind1(i,j))./area(i,j) - (dt.*dupwind1(i,j))./area(i,j);
                T2(i,j) = T2(i,j) + (dt.*upwind2(i,j))./area(i,j) - (dt.*dupwind2(i,j))./area(i,j);
                T3(i,j) = T3(i,j) + (dt.*upwind3(i,j))./area(i,j) - (dt.*dupwind3(i,j))./area(i,j);
                T4(i,j) = T4(i,j) + (dt.*upwind4(i,j))./area(i,j) - (dt.*dupwind4(i,j))./area(i,j);
                T5(i,j) = T5(i,j) + (dt.*upwind5(i,j))./area(i,j) - (dt.*dupwind5(i,j))./area(i,j);
                T6(i,j) = T6(i,j) + (dt.*upwind6(i,j))./area(i,j) - (dt.*dupwind6(i,j))./area(i,j);
                T7(i,j) = T7(i,j) + (dt.*upwind7(i,j))./area(i,j) - (dt.*dupwind7(i,j))./area(i,j);
                T8(i,j) = T8(i,j) + (dt.*upwind8(i,j))./area(i,j) - (dt.*dupwind8(i,j))./area(i,j);
            end
        end
    end

    B1 = T1(id);
    B2 = T2(id);
    B3 = T3(id);
    B4 = T4(id);
    B5 = T5(id);
    B6 = T6(id);
    B7 = T7(id);
    B8 = T8(id);
    
end