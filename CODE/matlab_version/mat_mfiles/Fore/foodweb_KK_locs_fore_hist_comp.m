% Create digraph objects for the food webs
% Forecast 

clear all
close all

datap = '/Volumes/FEISTY/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

%dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);

%% Hindcast
sname = 'Historic_';
load([dpath sname harv '_red_locs_biom_flux_KKfoodweb.mat']);

% Fluxes (consumption rates)
%rows 1:Npp->zoop, 2:zoop->F, 3:F->P, 4:NPP->det, 5:B->D, 6:P->D, 7:F->D 
%cols example, EBS, PUp, HOT
HedgesN(:,1) = [5e-1;5e-2;5e-3;5e-1;5e-2;5e-4;5e-3];
HedgesN(:,2) = flux(2,:)';
HedgesN(:,3) = flux(7,:)';
HedgesN(:,4) = flux(4,:)';
HedgesNC = num2cell(HedgesN);

HedgesC = {...
    'purple' 'yellow'   
    'yellow' 'red'     
    'red'    'blue'   
    'purple' 'brown'  
    'brown'  'green'  
    'blue'   'green'  
    'red'    'green'};

Hedges = [HedgesC HedgesNC];

% Biomasses
%rows 1:Npp, 2:zoop, 3:F, 4:P, 5:B, 6:D 
%cols example, EBS, PUp, HOT
Hxybio = [...
1 2  5   % 'purple'
2 2  10  % 'grey'  
3 2  2.5 % 'orange'
4 2  1   % 'blue'  
2 1  0.5 % 'black' 
4 1  1]; % 'green'

Hxybio(:,4) = bios(2,:)';
Hxybio(:,5) = bios(7,:)';
Hxybio(:,6) = bios(4,:)';

HG = cell(4,1);
for ii = 1:4
    HG{ii} = digraph(Hedges(:,1), Hedges(:,2), cell2mat(Hedges(:,ii+2)));
    HG{ii}.Nodes.bio = Hxybio(:,ii+2);
    HG{ii}.Nodes.x = Hxybio(:,1);
    HG{ii}.Nodes.y = Hxybio(:,2);
end


%% Forecast
sname = 'Forecast_';
load([dpath sname harv '_red_locs_biom_flux_KKfoodweb.mat']);

% Fluxes (consumption rates)
%rows 1:Npp->zoop, 2:zoop->F, 3:F->P, 4:NPP->det, 5:B->D, 6:P->D, 7:F->D 
%cols example, EBS, PUp, HOT
FedgesN(:,1) = [5e-1;5e-2;5e-3;5e-1;5e-2;5e-4;5e-3];
FedgesN(:,2) = flux(2,:)';
FedgesN(:,3) = flux(7,:)';
FedgesN(:,4) = flux(4,:)';
FedgesNC = num2cell(FedgesN);

FedgesC = {...
    'purple' 'yellow'   
    'yellow' 'red'     
    'red'    'blue'   
    'purple' 'brown'  
    'brown'  'green'  
    'blue'   'green'  
    'red'    'green'};

Fedges = [FedgesC FedgesNC];

% Biomasses
%rows 1:Npp, 2:zoop, 3:F, 4:P, 5:B, 6:D 
%cols example, EBS, PUp, HOT
Fxybio = [...
1 2  5   % 'purple'
2 2  10  % 'grey'  
3 2  2.5 % 'orange'
4 2  1   % 'blue'  
2 1  0.5 % 'black' 
4 1  1]; % 'green'

Fxybio(:,4) = bios(2,:)';
Fxybio(:,5) = bios(7,:)';
Fxybio(:,6) = bios(4,:)';

FG = cell(4,1);
for ii = 1:4
    FG{ii} = digraph(Fedges(:,1), Fedges(:,2), cell2mat(Fedges(:,ii+2)));
    FG{ii}.Nodes.bio = Fxybio(:,ii+2);
    FG{ii}.Nodes.x = Fxybio(:,1);
    FG{ii}.Nodes.y = Fxybio(:,2);
end

%% Combined
G(1) = HG(2);
G(2) = FG(2);
G(3) = HG(3);
G(4) = FG(3);
G(5) = HG(4);
G(6) = FG(4);

%% Some scaling setup

rmax = 0.5; % maximum radius 
bmax = 30;  % biomass represented by max radius
fmax = 1.0; % flux represented by max line width
wmax = 100;  % edge width scale

for ii = 1:6
    G{ii}.Nodes.r = sqrt((G{ii}.Nodes.bio .* (pi*rmax.^2)/30)/pi);
end

%--------------------
% Plot
%--------------------

th = linspace(0,2*pi,50);
nedge = numedges(G{1});
nnode = numnodes(G{1});

lbl = {'Hist EBS','Fore EBS','PUP','PUP','HOT','HOT'};

for ii = 1:6
    h(ii).ax = subplot(3,2,ii);

    % Plot edges first

    [h(ii).edg, Data] = plotdeb(G{ii}, 'p', 1/3, 'gmax', fmax, 'w', wmax);

    % Plot nodes

    h(ii).nd = arrayfun(@(x,y,r) patch(r.*cos(th)+x, r.*sin(th)+y, 'w'), ...
        G{ii}.Nodes.x, G{ii}.Nodes.y, G{ii}.Nodes.r, 'uni', 0);
    h(ii).nd = cat(1, h(ii).nd{:});

    % Color nodes, and match edge colors to their source nodes

    nvert = size(h(ii).edg.CData,1);

    set(h(ii).nd, 'facecolor', 'flat', 'edgecolor', 'w');
    set(h(ii).nd, {'CData'}, num2cell((1:nnode)'));
    h(ii).edg.CData = ones(nvert,1)*findnode(G{1}, G{1}.Edges.EndNodes(:,1))';

    title(lbl{ii});
end
set([h.ax], 'dataaspectratio', [1 1 1], ...
    'ylim', [0.5 2.75], 'xlim', [0.5 4.5], 'xtick', [], 'ytick', []);

%% Colors
% cmap = [...
%      0.49412      0.11765      0.61176   %purple
%      0.57255      0.58431      0.56863   %grey
%      0.97647      0.45098     0.023529   %orange
%      0.17255      0.43529      0.73333   %blue
%            0            0            0   %black
%     0.082353       0.6902      0.10196]; %green
cmap = [...
     0.57255      0.58431      0.56863   %grey
           1       0.8431            0   %yellow
     0.97647         0.19            0   %red
           0            0         0.65   %blue 
         0.4          0.2            0   %brown
         0.1         0.65      0.10196]; %green 
colormap(cmap);

set(gcf, 'color', 'w');
stamp(cfile)
print(gcf, '-dpng',[fpath 'Fore_Hist_foodweb_cbrt']);
