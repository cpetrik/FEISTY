%%% Update biomass
function [bio_sm,bio_md,predS,predM] = sub_update_be2_boris(BE,det,bio_in,Dcon,Dbio)
    %BE = benthic efficiency
    %det = seafloor detritus flux in g/m2/d
    %bio_in = benthic biomass in g/m2
    %Dcon = biomass specific consumption rate by MD & LD in /m2/d
    %Dbio = biomass of MD & LD in g/m2
    M = [0.02 11.2]; %mass in g
    
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
%     pars(1)=0.21;    % fraction of ingested material that is assimilated
%     pars(2)=0.61;    % respiration pre-factor
%     pars(3)=0.46e-2; % respiration scaling
%     pars(4)=1.76;    % ingestion pre-factor
%     pars(5)=-0.13;   % ingestion scaling
%     pars(7)=9.0e-4;  % mortality pre-factor
%     pars(8)=-0.40;   % mortality scaling
    
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
    
    bio_sm = Biom_out(:,1);
    bio_md = Biom_out(:,2);
    
    predS = eaten(:,1);
    predM = eaten(:,2);

end
