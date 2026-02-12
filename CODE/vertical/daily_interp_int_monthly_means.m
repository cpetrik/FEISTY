% Function for
% DAILY INTERPOLATION, INCLUDES UNIT CONVERSION

function ESM = daily_interp_int_monthly_means(NWID,Time,Tdays,...
    TP,TB,Det,MdZ,LgZ,MZl,LZl,thk,WID,ni,nj)

%% Integrate 3-d first
TEMP_z = squeeze(sum(TP .* thk,3,'omitnan') ./ sum(thk,3,'omitnan'));
MZ = squeeze(sum(MdZ .* thk,3,'omitnan'));
LZ = squeeze(sum(LgZ .* thk,3,'omitnan'));
MZloss = squeeze(sum(MZl .* thk,3,'omitnan'));
LZloss = squeeze(sum(LZl .* thk,3,'omitnan'));

%% Interpolate
ESM.Tp = nan(NWID,length(Tdays));
ESM.Zm = nan(NWID,length(Tdays));
ESM.Zl = nan(NWID,length(Tdays));
ESM.dZm = nan(NWID,length(Tdays));
ESM.dZl = nan(NWID,length(Tdays));
ESM.dz = nan(NWID,length(Tdays));

for W = 1:NWID

    [m,n] = ind2sub([ni,nj],WID(W)); % spatial index of water cell
    % location of interest
    Tp  = double(squeeze(TEMP_z(m,n,:)));
    Zm  = double(squeeze(MZ(m,n,:)));
    Zl  = double(squeeze(LZ(m,n,:)));
    dZm = double(squeeze(MZloss(m,n,:)));
    dZl = double(squeeze(LZloss(m,n,:)));
    Tb  = double(squeeze(TB(m,n,:)));
    det = double(squeeze(Det(m,n,:)));

    Tb = Tb';
    det = det';

    % bottom temperature (in Celcius)
    tb = interp1(Time, Tb, Tdays,'linear','extrap');
    ESM.Tb(W,:) = tb;

    % detrital flux to benthos: 'molN m-2 s-1' to g(WW) m-2 d-1
    % 106/16 C to N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    % 60s * 60m * 24h
    de = interp1(Time, det, Tdays,'linear','extrap');
    ESM.det(W,:) = de * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

    % pelagic temperature (in Celcius)
    tp = interp1(Time, Tp, Tdays,'linear','extrap');
    ESM.Tp(W,:) = tp;

    % medium zoo: 'molN/kg' to g(WW) m-3
    % 1 m3 = 1035 kg
    % 106/16 C to N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    mz = interp1(Time, Zm, Tdays,'linear','extrap');
    ESM.Zm(W,:) = mz * 1035 * (106/16) * 12.01 * 9.0;

    % large zoo: 'molN/kg' to g(WW) m-3
    lz = interp1(Time, Zl, Tdays,'linear','extrap');
    ESM.Zl(W,:) = lz * 1035 * (106/16) * 12.01 * 9.0;

    % med zoo mortality: 'mol N kg-1 s-1' to g(WW) m-3 d-1
    % 1 m3 = 1035 kg
    % 106/16 C to N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W
    % 60s * 60m * 24h
    mzl = interp1(Time, dZm, Tdays,'linear','extrap');
    ESM.dZm(W,:) = mzl * 1035 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;

    % lrg zoo mortality: 'mol N kg-1 s-1' to g(WW) m-3 d-1
    lzl = interp1(Time, dZl, Tdays,'linear','extrap');
    ESM.dZl(W,:) = lzl * 1035 * (106/16) * 12.01 * 9.0 * 60 * 60 * 24;


end

% Negative biomass or mortality loss from interp
ESM.Zm(ESM.Zm<0) = 0.0;
ESM.Zl(ESM.Zl<0) = 0.0;
ESM.dZm(ESM.dZm<0) = 0.0;
ESM.dZl(ESM.dZl<0) = 0.0;
ESM.det(ESM.det<0) = 0.0;

ESM.det(isnan(ESM.det)) = 0.0;
ESM.Zm(isnan(ESM.Zm)) = 0.0;
ESM.Zl(isnan(ESM.Zl)) = 0.0;
ESM.dZm(isnan(ESM.dZm)) = 0.0;
ESM.dZl(isnan(ESM.dZl)) = 0.0;


end
