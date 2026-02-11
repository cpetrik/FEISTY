% Function for
% DAILY INTERPOLATION, INCLUDES UNIT CONVERSION

function ESM = daily_interp_monthly_means(nz,Time,Tdays,...
    Tp,Tb,det,Zm,Zl,dZm,dZl,thk)


%btm vars

% bottom temperature (in Celcius)
tb = interp1(Time, Tb, Tdays,'linear','extrap');
ESM.Tb = tb;

% detrital flux to benthos: 'molN m-2 s-1' to g(WW) m-2 d-1 
% 106/16 C to N
% 12.01 g C in 1 mol C
% 1 g dry W in 9 g wet W (Pauly & Christiansen)
% 60s * 60m * 24h
de = interp1(Time, det, Tdays,'linear','extrap');
ESM.det = de * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;
ESM.det(ESM.det<0) = 0.0;
ESM.det(isnan(ESM.det)) = 0.0;

ESM.Tp = nan(nz,length(Tdays));
ESM.Zm = nan(nz,length(Tdays));
ESM.Zl = nan(nz,length(Tdays));
ESM.dZm = nan(nz,length(Tdays));
ESM.dZl = nan(nz,length(Tdays));
ESM.dz = nan(nz,length(Tdays));

for z=1:nz
    % pelagic temperature (in Celcius)
    tp = interp1(Time, Tp(z,:), Tdays,'linear','extrap');
    ESM.Tp(z,:) = tp;

    % medium zoo: 'molN/kg' to g(WW) m-3
    % 1 m3 = 1035 kg
    % 106/16 C to N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    mz = interp1(Time, Zm(z,:), Tdays,'linear','extrap');
    ESM.Zm(z,:) = mz * 1035 * (106/16) * 12.01 * 9.0;

    % large zoo: 'molN/kg' to g(WW) m-3
    lz = interp1(Time, Zl(z,:), Tdays,'linear','extrap');
    ESM.Zl(z,:) = lz * 1035 * (106/16) * 12.01 * 9.0;

    % med zoo mortality: 'mol N kg-1 s-1' to g(WW) m-3 d-1
    % 1 m3 = 1035 kg
    % 106/16 C to N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W 
    % 60s * 60m * 24h
    mzl = interp1(Time, dZm(z,:), Tdays,'linear','extrap');
    ESM.dZm(z,:) = mzl * 1035 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

    % lrg zoo mortality: 'mol N kg-1 s-1' to g(WW) m-3 d-1
    lzl = interp1(Time, dZl(z,:), Tdays,'linear','extrap');
    ESM.dZl(z,:) = lzl * 1035 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

    % grid cell thickness
    tk = interp1(Time, thk(z,:), Tdays,'linear','extrap');
    ESM.dz(z,:) = tk;

end

% Negative biomass or mortality loss from interp
ESM.Zm(ESM.Zm<0) = 0.0;
ESM.Zl(ESM.Zl<0) = 0.0;
ESM.dZm(ESM.dZm<0) = 0.0;
ESM.dZl(ESM.dZl<0) = 0.0;

ESM.Zm(isnan(ESM.Zm)) = 0.0;
ESM.Zl(isnan(ESM.Zl)) = 0.0;
ESM.dZm(isnan(ESM.dZm)) = 0.0;
ESM.dZl(isnan(ESM.dZl)) = 0.0;

ESM.dz(ESM.dz<0) = nan;


end
