%% BORIS benthic biomass eqs
clear all close all

% Set up variables & params
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');
det=mean(COBALT.det,2);

load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;
%%
%! Fish biomass
X = 1.0;
Md.bio = ones(NX,1) * X;
Ld.bio = ones(NX,1) * X;

%! Rates, etc.
nzero = {'con_be' 'I' 'pred'};
for i = 1 : length(nzero)
    Md.(nzero{i}) = zeros(NX,1);
    Ld.(nzero{i}) = zeros(NX,1);
end

%! Detritus
BENT.sm = ones(NX,1) * X;
BENT.md = ones(NX,1) * X;
BENT.predS = zeros(NX,1);
BENT.predM = zeros(NX,1);

% BORIS orig masses in gWW
%m = [8.9088E-07	1.78176E-06	3.56352E-06	7.12704E-06	1.42541E-05	2.85082E-05	5.70163E-05	0.000114033	0.000228065	0.000456131	0.000912261	0.001824522	0.003649044	0.007298089	0.014596178	0.029192356];
%mlog = log10(m); dmlog = diff(mlog); %=0.3010
% FEISTY Small mean mass = 0.02 g; Medium mean mass = 11.2 g
%fish = [0.02 11.2]; flog = log10(fish); dflog = diff(flog); %=2.7482
M = [0.02 11.2];

%BE = assim * (1 - fbact) = 0.1 * (1-0.25)
%but explicitly model assim below, so remove it here
BE = 0.75;

bio_in = [BENT.sm,BENT.md];
Dcon = [Md.con_be,Ld.con_be] + 1e-5;
Dbio = [Md.bio,Ld.bio];

%% Update biomass
%function [bio_sm,bio_md,predS,predM] = sub_update_be2_boris(BE,det,bio_in,Dcon,Dbio)
%BE = benthic efficiency
%det = seafloor detritus flux in g/m2/d
%bio_in = benthic biomass in g/m2
%Dcon = biomass specific consumption rate by MD & LD in /m2/d
%Dbio = biomass of MD & LD in g/m2

% gi        = Specific ingestion rate (d-1) for each body size class
% mi        = Specific mortality rate (d-1)
% ri        = Specific respiration rate (d-1)
% kingi     = Ingestion rate half saturation constant (g R m-2)
% kmi       = Mortality rate half saturation constant (g B m-2)
% alf_i     = assimilation efficiency
% r         = Input of detritus (POC flux) (g m-2 d-1)

%! Universal optimal parameters from Andrew Y.
pars(1)=0.4;    % fraction of ingested material that is assimilated
pars(2)=0.53;   % respiration pre-factor
pars(3)=0.03;   % respiration scaling
pars(4)=0.02;   % ingestion pre-factor
pars(5)=-0.1;   % ingestion scaling
pars(7)=0.0017; % mortality pre-factor
pars(8)=-0.37;  % mortality scaling

%! Optimized parameters in Yool et al. 2017
pars(1)=0.21;    % fraction of ingested material that is assimilated
pars(2)=0.61;    % respiration pre-factor
pars(3)=0.46e-2; % respiration scaling
pars(4)=1.76;    % ingestion pre-factor
pars(5)=-0.13;   % ingestion scaling
pars(7)=9.0e-4;  % mortality pre-factor
pars(8)=-0.40;   % mortality scaling

%For mass conservation at steady state, total respiration must equal
%the flux of detritus from the overlying water i.e.
% sum(ri*(alf_i*G_i)*Bi) = (BE)*r

%! Mortality from predation
eaten = Dcon .* Dbio;

%! Detritus available to benthos
r = BE .* det; %Needs to be in units of per time (g/m2/d) * (g/m2)

%! Physiology params
% assimilation efficiencies for detritus
assimcff = pars(1);
% respiration coefficient (range 0.1-0.9)
resp_cff = pars(2) * M.^pars(3);
% specific ingestion rate  (/day)
gi  = pars(4) * M.^pars(5);
% specific mortality rate  (/day)
mi = pars(7) * M.^pars(8);
%! Physiology eqs
I        = gi;
ingest   = I .* bio_in .* r;           % ingestion rate
assim    = assimcff .* ingest;         % assimilation rate
%defec    = (1-assimcff) .* ingest;    % defecation rate
resp     = resp_cff .* assim;          % respiration rate
%death    =  mi.* bio_in.*bio_in;      % density dependent mortality rate
net_prod = assim-resp;

death = eaten;

Biom_out =  assim - resp - death;      % biomass equation

%R_out   =  source + sum(defec) + sum(death)- sum(ingest); % detritus equation

predS = eaten(:,1);
predM = eaten(:,2);

%% Plot physiol in comp to fish
temp = 10;
cmax = (exp(0.063*(temp-10.0)) .* 20 .* M.^(-0.25)) ./365.0;
met = 0.2 * (exp(0.0905*(temp-10.0)) .* 20 .* M.^(-0.175)) ./365.0;
enc = (exp(0.063*(temp-10.0)) .* 70 .* M.^(-0.2)) ./365.0;

%figure(1)
figure
plot(log10(M),log10(gi),'.-r','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(M),log10(cmax),'.-b','MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 ingestion (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Bent','Fish')
legend('location','southwest')
title('Specific ingestion rate')

%figure(2)
figure
plot(log10(M),log10(resp_cff),'.-r','MarkerSize',15,'LineWidth',2); hold on;
plot(log10(M),log10(met),'.-b','MarkerSize',15,'LineWidth',2); hold on;
ylabel('log10 respiration (gWW gWW^-^1 d^-^1)')
xlabel('log10 weight (gWW)')
legend('Bent','Fish')
legend('location','southwest')
title('Specific respiration rate')


