% Standalone code for running BORIS

clear all
close all

%% Model time and place
%! How long to run the model
YEARS = 1;%50;
DAYS = 365*YEARS;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
DT = 1;

%! Grid; choose where and when to run the model
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;
% Dims
nt = 12;

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');


%% Storage variables
S_Bent_Sbio = zeros(NX,DAYS);
S_Bent_Mbio = zeros(NX,DAYS);
S_Bent_Smort = zeros(NX,DAYS);
S_Bent_Mmort = zeros(NX,DAYS);
S_Det = zeros(NX,DAYS);

Spinup_SBent.bio = NaN*ones(NX,nt);
Spinup_MBent.bio = NaN*ones(NX,nt);
Spinup_SBent.mort = NaN*ones(NX,nt);
Spinup_MBent.mort = NaN*ones(NX,nt);
Spinup_Det = NaN*ones(NX,nt);

X = 1e-2;
BENT.sm = ones(NX,1) * X;
BENT.md = ones(NX,1) * X;
BENT.mortS = zeros(NX,1);
BENT.mortM = zeros(NX,1);
Det = zeros(NX,1);

%% Parameters
% Transfer of detritus to benthos (not bacteria)
bent_eff = 0.75;

% Background mortality
Nat_mrt = 0.1/365;

%mass in g
M = [0.02 11.2];

% gi        = Specific ingestion rate (d-1) for each body size class
% mi        = Specific mortality rate (d-1)
% ri        = Specific respiration rate (d-1)
% kingi     = Ingestion rate half saturation constant (g R m-2)
% kmi       = Mortality rate half saturation constant (g B m-2)
% alf_i     = assimilation efficiency
% r         = Input of detritus (POC flux) (g m-2 d-1)

%! Universal optimal parameters from Andrew Y.
% pars(1)=0.4;    % fraction of ingested material that is assimilated
% pars(2)=0.53;   % respiration pre-factor
% pars(3)=0.03;   % respiration scaling
% pars(4)=0.02;   % ingestion pre-factor
% pars(5)=-0.1;   % ingestion scaling
% pars(7)=0.0017; % mortality pre-factor
% pars(8)=-0.37;  % mortality scaling
%with these params, resp rate is 0.1885 g/m2/d  &  0.2279 g/m2/d for my
%sizes
%but mean of COBALT bottom detritus is << that  
%quantile(COBALT.det(:),[0.05 0.25 0.5 0.75 0.95])
%               0.0012  0.0074  0.0175  0.0402  0.8322

%! Optimized parameters in Yool et al. 2017
pars(1)=0.21;    % fraction of ingested material that is assimilated
pars(2)=0.61;    % respiration pre-factor
pars(3)=0.46e-2; % respiration scaling
pars(4)=1.76;    % ingestion pre-factor
pars(5)=-0.13;   % ingestion scaling
pars(7)=9.0e-4;  % mortality pre-factor
pars(8)=-0.40;   % mortality scaling
%with these params, resp rate is 0.1258 g/m2/d  &  0.1295 g/m2/d for my
%sizes
%on day 1, 8% and 5% had positive gain, all others died

%! Physiology params
%(What if rates are per yr instead of per day?) - NO don't eat enough
% assimilation efficiencies for detritus
assimcff = pars(1);
% respiration coefficient (range 0.1-0.9)
resp_cff = pars(2) * M.^pars(3);
% specific ingestion rate  (/day)
gi  = pars(4) * M.^pars(5);
% specific mortality rate  (/day)
mi = pars(7) * M.^pars(8);


%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Run model with no fishing
MNT=0;
for YR = 1:YEARS % years
    num2str(YR)
    
    for DAY = 1:20;%:DT:DAYS % days
        
        % Future time step
        DY = int64(ceil(DAY));
        
        % Call COBALT detritus
        ENVR = get_COBALT(COBALT,ID,DY);
        ENVR.det = sub_neg(ENVR.det);
        
        %% Benthos
        %For mass conservation at steady state, total respiration must equal
        %the flux of detritus from the overlying water i.e.
        % sum(ri*(alf_i*G_i)*Bi) = (BE)*r
        
        %! Mortality from predation
        %eaten = Dcon .* Dbio;
        eaten = 0;
        
        %! Detritus available to benthos
        adet = bent_eff .* ENVR.det; %Needs to be in units of per time (g/m2/d) * (g/m2)
        r = Det + adet;
        
        %! Physiology eqs 
        bio_in = [BENT.sm BENT.md];
        I        = gi;                         % ingestion rate
        ingest   = I .* bio_in .* r;           % ingested biomass
        assim    = assimcff .* ingest;         % assimilated biomass
        defec    = (1-assimcff) .* ingest;     % defecation rate
        %resp     = resp_cff .* assim;          % respiration rate
        resp     = resp_cff .* assimcff .* bio_in; % resp rate should NOT depend on amount assimilated
        %death    =  mi.* bio_in.*bio_in;      % density dependent mortality rate
        net_prod = assim - resp;
        
        %death = Nat_mrt.* bio_in + eaten;
        death = Nat_mrt.* bio_in;
        
        Biom_out =  assim - resp - death;      % biomass equation
        
        % What if you need to accumulate det?
        %R_out   =  source + sum(defec) + sum(death)- sum(ingest); % detritus equation
        %Det = adet + sum(defec,2) + sum(death,2) - sum(ingest,2);
        
        
        BENT.sm = Biom_out(:,1);
        BENT.md = Biom_out(:,2);
        
        BENT.sm = sub_check(BENT.sm);
        BENT.md = sub_check(BENT.md);
        %Det = sub_check(Det);
        
        BENT.mortS(ID) = net_prod(:,1);
        BENT.mortM(ID) = net_prod(:,2);
        
        %% Save last year
        %if (YR==YEARS)
        S_Bent_Sbio(:,DY) = BENT.sm;
        S_Bent_Mbio(:,DY) = BENT.md;
        
        S_Bent_Smort(:,DY) = BENT.mortS;
        S_Bent_Mmort(:,DY) = BENT.mortM;
        
        %end
        
    end %Days
    
end %Years

%%
figure
subplot(2,1,1)
plot(1:20,mean(S_Bent_Sbio(:,1:20)),'b'); hold on;

subplot(2,1,2)
plot(1:20,mean(S_Bent_Mbio(:,1:20)),'k')

%one location seems to exponentially incr and throws off the mean

%% Calculate monthly means and save
aa = (cumsum(MNTH)+1);
a = [1,aa(1:end-1)]; % start of the month
b = cumsum(MNTH); % end of the month
for i = 1:12
    %! Put vars of netcdf file
    Spinup_SBent.bio(:,i) = mean(S_Bent_Sbio(:,a(i):b(i)),2);
    Spinup_MBent.bio(:,i) = mean(S_Bent_Mbio(:,a(i):b(i)),2);
    
    Spinup_SBent.mort(:,i) = mean(S_Bent_Smort(:,a(i):b(i)),2);
    Spinup_MBent.mort(:,i) = mean(S_Bent_Mmort(:,a(i):b(i)),2);
    
end %Monthly mean

%% Save
% save([fname,'.mat'],...
%     'Spinup_Sml_d','Spinup_Med_d','Spinup_Lrg_d',...
%     'Spinup_SBent','Spinup_MBent');






