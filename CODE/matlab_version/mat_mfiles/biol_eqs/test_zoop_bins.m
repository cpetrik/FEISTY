%Example of dividing zooplankton into multiple size bins

%Mesozoop biom ranges from 1e2 to 1e4 mg C/m2

%BMZ = 3.6; 5.8 Pg C/yr =1e15gC/yr
%BLZ = 1.1; 1.7 Pg C/yr

%%
%365d/yr
%360e6 km2
%perkm2 = 3.6e15 / 365 / 360e6;
%perm2 = perkm2 / 1e6;

BMZ = 3.6e15 / 365 / 360e6; %/ 1e6; 
BLZ = 1.1e15 / 365 / 360e6; %/ 1e6;

%%

mlow = 0.2;
mhi = 2;
llow = 2; 
lhi = 20;
hm = 0.2;
hl = 2;
mbin_ends = mlow:hm:mhi;
mbin_mids = ((mbin_ends(1)+mbin_ends(2))/2):hm:((mbin_ends(end-1)+mbin_ends(end))/2);
lbin_ends = llow:hl:lhi;
lbin_mids = ((lbin_ends(1)+lbin_ends(2))/2):hl:((lbin_ends(end-1)+lbin_ends(end))/2);
Hm = mbin_ends(end) - mbin_ends(1);
Hl = lbin_ends(end) - lbin_ends(1);
Mmed = (mlow+mhi)/2;
Lmed = (llow+lhi)/2;

%% Calc biom/mm slope
s = (BLZ-BMZ) / (Lmed-Mmed);

%% Solve for a1
aM1 = (BMZ - 0.5*Hm^2 * s) / Hm;
aL1 = (BLZ - 0.5*Hl^2 * s) / Hl;

%% Calc biomass in each bin
BM = NaN*ones(length(mbin_mids),1);
BL = NaN*ones(length(lbin_mids),1);
BM(1) = aM1*hm + 0.5*hm^2 * s;
BL(1) = aL1*hl + 0.5*hl^2 * s;

for i=2:length(lbin_mids)
    BM(i) = aM1*hm + 0.5*hm^2 * s + (i-1)*hm^2 * s;
    BL(i) = aL1*hl + 0.5*hl^2 * s + (i-1)*hl^2 * s;
end

%% If a1 is neg, move forward until find pos
j=0:(length(lbin_mids)-1);
aL = aL1 + j.*hl.*s;
pid = find(aL>0);
start=pid(1)-1;

Hl = hl*length(start:length(lbin_mids));
aL2 = NaN*ones(length(lbin_mids),1);
aL2(start) = (BLZ - 0.5*Hl^2 * s) / Hl;


BL2 = zeros(length(lbin_mids),1);
BL2(start) = aL2(start)*hl + 0.5*hl^2 * s;
for i=(start+1):length(lbin_mids)
    BL2(i) = BL2(start)*hl + 0.5*hl^2 * s + (i-1)*hl^2 * s;
end

% THIS DOES NOT WORK


