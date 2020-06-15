%%%%% Fishing
function [MFbio, LPbio, LDbio] = sub_fishing_mass(MFbio,LPbio,LDbio,AREA)
	if FISHING > 0.0
        ALL_pl = MFbio .* AREA;
		ALL_pi = LPbio .* AREA;
		ALL_de = LDbio .* AREA;

		%! Total fish biomass
		TOT = sum(ALL_pi) + sum(ALL_pl) + sum(ALL_de);
        ALL_pl = (ALL_pl./TOT) .* FISHING;
		ALL_pi = (ALL_pi./TOT) .* FISHING;
		ALL_de = (ALL_de./TOT) .* FISHING;

		%! Calc total biomass of fish in the ocean
        MFbio = ALL_pl ./ AREA;
		LPbio = ALL_pi ./ AREA;
		LDbio = ALL_de ./ AREA;
	end
	%return MFbio, LPbio, LDbio
end