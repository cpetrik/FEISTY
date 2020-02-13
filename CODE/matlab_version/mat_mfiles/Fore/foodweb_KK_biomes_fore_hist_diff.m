% Create digraph objects for the food webs
% Forecast 

clear all
close all

datap = '/Volumes/FEISTY/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

%dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' dp '/'];
lpath=['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Data/' dp '/'];

ppath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);

%% Hindcast
load([fpath 'Hist_All_fish03_biomes_biom_flux_KKfoodweb_prods.mat']);

% Fluxes (consumption rates)
%rows 1:Npp->zoop, 2:zoop->F, 3:F->P, 4:NPP->det, 5:B->D, 6:P->D, 7:F->D 
%cols example, EBS, PUp, HOT
HedgesN(:,1) = [5e-1;5e-2;5e-3;5e-1;5e-2;5e-4;5e-3];
HedgesN(:,2) = flux(1,:)';
HedgesN(:,3) = flux(2,:)';
HedgesN(:,4) = flux(3,:)';
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

Hxybio(:,4) = bios(1,:)';
Hxybio(:,5) = bios(2,:)';
Hxybio(:,6) = bios(3,:)';

%% Forecast
load([fpath 'Fore_All_fish03_biomes_biom_flux_KKfoodweb_prods.mat']);

% Fluxes (consumption rates)
%rows 1:Npp->zoop, 2:zoop->F, 3:F->P, 4:NPP->det, 5:B->D, 6:P->D, 7:F->D 
%cols example, EBS, PUp, HOT
FedgesN(:,1) = [5e-1;5e-2;5e-3;5e-1;5e-2;5e-4;5e-3];
FedgesN(:,2) = flux(1,:)';
FedgesN(:,3) = flux(2,:)';
FedgesN(:,4) = flux(3,:)';
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

Fxybio(:,4) = bios(1,:)';
Fxybio(:,5) = bios(2,:)';
Fxybio(:,6) = bios(3,:)';

%% Combined
edgesN = HedgesN - FedgesN;
xybio(:,1:2) = Hxybio(:,1:2);
xybio(:,3:5) = Hxybio(:,4:6) - Fxybio(:,4:6);

% Just decr
neg1 = find(edgesN(:)<0);
edgesN(neg1) = 0;
neg2 = find(xybio(:)<0);
xybio(neg2) = 0;

edgesNC = num2cell(edgesN);
edgesC = {...
    'purple' 'yellow'   
    'yellow' 'red'     
    'red'    'blue'   
    'purple' 'brown'  
    'brown'  'green'  
    'blue'   'green'  
    'red'    'green'};
edges = [edgesC edgesNC];

G = cell(3,1);
for ii = 1:3
    G{ii} = digraph(edges(:,1), edges(:,2), cell2mat(edges(:,ii+2)));
    G{ii}.Nodes.bio = xybio(:,ii+2);
    G{ii}.Nodes.x = xybio(:,1);
    G{ii}.Nodes.y = xybio(:,2);
end

%% Some scaling setup

rmax = 0.5; % maximum radius 
bmax = 30;  % biomass represented by max radius
fmax = 1.0; % flux represented by max line width
wmax = 100;  % edge width scale

for ii = 1:3
    G{ii}.Nodes.r = sqrt((G{ii}.Nodes.bio .* (pi*rmax.^2)/30)/pi);
end

%--------------------
% Plot
%--------------------

th = linspace(0,2*pi,50);
nedge = numedges(G{1});
nnode = numnodes(G{1});

lbl = bnames;

for ii = 1:3
    h(ii).ax = subplot(3,1,ii);

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

% Colors
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
print(gcf, '-dpng',[ppath 'Fore_Hist_foodweb_cbrt_biomes_diff_decr']);


%% Combined incr
edgesN = FedgesN - HedgesN;
xybio(:,1:2) = Hxybio(:,1:2);
xybio(:,3:5) = Fxybio(:,4:6) - Hxybio(:,4:6);

% Just incr
neg1 = find(edgesN(:)<0);
edgesN(neg1) = 0;
neg2 = find(xybio(:)<0);
xybio(neg2) = 0;

edgesNC = num2cell(edgesN);
edgesC = {...
    'purple' 'yellow'   
    'yellow' 'red'     
    'red'    'blue'   
    'purple' 'brown'  
    'brown'  'green'  
    'blue'   'green'  
    'red'    'green'};
edges = [edgesC edgesNC];

G = cell(3,1);
for ii = 1:3
    G{ii} = digraph(edges(:,1), edges(:,2), cell2mat(edges(:,ii+2)));
    G{ii}.Nodes.bio = xybio(:,ii+2);
    G{ii}.Nodes.x = xybio(:,1);
    G{ii}.Nodes.y = xybio(:,2);
end

% Some scaling setup

rmax = 0.5; % maximum radius 
bmax = 30;  % biomass represented by max radius
fmax = 1.0; % flux represented by max line width
wmax = 100;  % edge width scale

for ii = 1:3
    G{ii}.Nodes.r = sqrt((G{ii}.Nodes.bio .* (pi*rmax.^2)/30)/pi);
end

%--------------------
% Plot
%--------------------
close all

th = linspace(0,2*pi,50);
nedge = numedges(G{1});
nnode = numnodes(G{1});


for ii = 1:3
    h(ii).ax = subplot(3,1,ii);

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

% Colors
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
print(gcf, '-dpng',[ppath 'Fore_Hist_foodweb_cbrt_biomes_diff_incr']);
