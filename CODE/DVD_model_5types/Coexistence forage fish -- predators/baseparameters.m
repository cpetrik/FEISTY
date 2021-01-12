clearvars param y
param.tEnd =  300;
%
% Resources:
%
param.ixR = [1 2 3 4];
param.w(param.ixR) = [2e-06  0.001 0.5e-03 0.25]; % lower limit
param.wc(param.ixR) = [2e-06*sqrt(500) 0.001*sqrt(500) 0.25 0.25*sqrt(500)]; % central size
param.wu(param.ixR) = [0.001 0.5  125  125]; % upper limit
param.r =  [1  1  1  1];  % g ww/m2/yr
param.K =  [5  5  5  0];  % g ww/m2
%
% stages
%   
param.nstage = 3; % number of stages predator use 3, 6, 9, etc (prey = 2/3)
param.nsize  = param.nstage + 1; % 
param.sfish = 0.001; % smallest size fish (all fish)
param.lfish = 125000; % largest size fish (only predator)
param.smat = 0.5; % weight at maturity forage/meso
param.lmat = 250; % weight at maturity predators
param.sizes = logspace(log10(param.sfish), log10(param.lfish),param.nsize);
[nb,param.maxsmall] = min(abs(param.sizes-250));
param.stage = param.sizes(end)/param.sizes(end-1);

%
% Species:
%
param.nSpecies = 5;
param.ixFish = 5: 4+(param.nstage*2/3)*2+(param.nstage)*3; % Index for all fish
%
%
% Indices and weight classes small pelagics
param.w = [param.w param.sizes(1:param.maxsmall-1)];
param.ix1(1) = 4+1;
param.ix2(1) = 4+(param.maxsmall-1);
% Indices and weight classes mesopelagics:
param.w = [param.w param.sizes(1:param.maxsmall-1)];
param.ix1(2) = param.ix2(1) + 1; 
param.ix2(2) = param.ix2(1) + (param.maxsmall-1);
% Indices and weight classes large pelagic:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(3) = param.ix2(2) + 1; 
param.ix2(3) = param.ix2(2) + (param.nsize-1);
% Indices and weight classes bathypelagic:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(4) = param.ix2(3) + 1; 
param.ix2(4) = param.ix2(3) + (param.nsize-1);
% Indices and weight classes large demersal:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(5) = param.ix2(4) + 1; 
param.ix2(5) = param.ix2(4) + (param.nsize-1);

param.wc(param.ixFish) = param.w(param.ixFish)*sqrt(param.stage); % central sizes
param.wu(param.ixFish) = param.w(param.ixFish)*param.stage; % Upper sizes

% predator prey preference
param.beta = 400; 
param.sigma = 1.3;
%
% Initial conditions 
%
param.B0 = 0*param.ixFish+.1;
%param.B0(param.ix1(1)-4:param.ix2(1)-4)=0; %put small pelagics to zero
%param.B0(param.ix1(2)-4:param.ix2(2)-4)=0; %put small mesopelagics to zero
%param.B0(param.ix1(3)-4:param.ix2(3)-4)=0; %put large pelagics to zero
%param.B0(param.ix1(4)-4:param.ix2(4)-4)=0; %put bathypelagics to zero
%param.B0(param.ix1(5)-4:param.ix2(5)-4)=0; %put large demersals to zero

%
% Habitat and interactions
%
param.bottom = 100; % depth in meter
param.photic = 150;
param.mesop = 250;
param.visual = 1.5; % scalar; >1 visual predation primarily during the day, = 1 equal day and night
[param.theta, param.depthDay, param.depthNight, param.avlocDay, param.avlocNight] = calcpreference(param); % feeding preference matrix 

%
% Calculate temperature dependency 
%
param.Q10 = 1.88;
param.region = 4; % 1=tropical, 2=temperate, 3=boreal, 4 = no temp scaling
[param.scTemp] = calctemperature(param);

%
% Physiology:
%
param.h = 20; 
param.met = 0.2*param.h; % maintenance costs, 20% of param.h
param.epsAssim = 0.7;
param.q = 0.8;
param.n = 0.75;
param.m = 0.825;
param.F = 0*param.ixFish';  % Fishing mortality
param.gamma = 70; % factor for the max clearance rate (area per time) 
param.eRepro = repmat(0.01,param.nSpecies,1)';
param.mort0 = 0*param.ixFish'+.1; 

% define investment in maturation
[~,param.matstageS] = min(abs(param.sizes-param.smat));
[~,param.matstageL] = min(abs(param.sizes-param.lmat));
param.kappaS = [ones(param.matstageS-1,1)' repmat(0.5,(param.maxsmall-1-(param.matstageS-1)),1)'];
param.kappaL = [ones(param.matstageL-1,1)' repmat(0.5,(param.nsize-1-(param.matstageL-1)),1)'];
param.kappa = [param.kappaS param.kappaS param.kappaL param.kappaL param.kappaL];

param.z = (param.w./param.wu)';
param.Cmax = (param.h*param.wc.^param.n)./param.wc .* param.scTemp;
param.V = (param.gamma*param.wc.^param.q)./param.wc;
param.Mc = (param.met*param.wc.^param.m)./param.wc .* param.scTemp;