% Create digraph objects for the food webs

clear all
close all

datap = '/Volumes/GFDL/NC/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

%dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
dp = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
sname = 'Climatol_';
harv = 'All_fish03';
dpath = [datap char(dp) '/'];
fpath = [figp char(dp) '/'];
if (~isdir([figp char(dp)]))
    mkdir([figp char(dp)])
end
cfile = char(dp);
load([dpath sname harv '_red_locs_biom_flux_KKfoodweb.mat']);

% Fluxes (consumption rates)
%rows 1:Npp->zoop, 2:zoop->F, 3:F->P, 4:NPP->det, 5:B->D, 6:P->D, 7:F->D 
%cols example, EBS, PUp, HOT
edgesN(:,1) = [5e-1;5e-2;5e-3;5e-1;5e-2;5e-4;5e-3];
edgesN(:,2) = flux(2,:)';
edgesN(:,3) = flux(7,:)';
edgesN(:,4) = flux(4,:)';
edgesNC = num2cell(edgesN);

% edgesC = {...
%     'purple' 'grey'   
%     'grey'   'orange'     
%     'orange' 'blue'   
%     'purple' 'black'  
%     'black'  'green'  
%     'blue'   'green'  
%     'orange' 'green'};
edgesC = {...
    'purple' 'yellow'   
    'yellow' 'red'     
    'red'    'blue'   
    'purple' 'brown'  
    'brown'  'green'  
    'blue'   'green'  
    'red'    'green'};

edges = [edgesC edgesNC];

% Biomasses
%rows 1:Npp, 2:zoop, 3:F, 4:P, 5:B, 6:D 
%cols example, EBS, PUp, HOT
xybio = [...
1 2  5   % 'purple'
2 2  10  % 'grey'  
3 2  2.5 % 'orange'
4 2  1   % 'blue'  
2 1  0.5 % 'black' 
4 1  1]; % 'green'

xybio(:,4) = bios(2,:)';
xybio(:,5) = bios(7,:)';
xybio(:,6) = bios(4,:)';

G = cell(4,1);
for ii = 1:4
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

for ii = 1:4
    G{ii}.Nodes.r = sqrt((G{ii}.Nodes.bio .* (pi*rmax.^2)/30)/pi);
end

%--------------------
% Plot
%--------------------

th = linspace(0,2*pi,50);
nedge = numedges(G{1});
nnode = numnodes(G{1});

lbl = {'Reference sizes','EBS', 'PUP', 'HOT'};

for ii = 1:4
    h(ii).ax = subplot(2,2,ii);

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
%
subplot(2,2,1)
%Numbers
text(1,2.6,sprintf('%2.1f',xybio(1,3)),'HorizontalAlignment','center'); hold on;
text(2,2.6,sprintf('%2.1f',xybio(2,3)),'HorizontalAlignment','center'); hold on;
text(2,0.6,sprintf('%2.1f',xybio(5,3)),'HorizontalAlignment','center'); hold on;
text(3,2.6,sprintf('%2.1f',xybio(3,3)),'HorizontalAlignment','center'); hold on;
text(4,2.6,sprintf('%2.1f',xybio(4,3)),'HorizontalAlignment','center'); hold on;
text(4,0.6,sprintf('%2.1f',xybio(6,3)),'HorizontalAlignment','center'); hold on;
text(1.5,2.3,sprintf('%2.1f',edges{1,3}),'HorizontalAlignment','center'); hold on;
text(2.5,2.3,sprintf('%1.0e',edges{2,3}),'HorizontalAlignment','center'); hold on;
text(3.5,2.3,sprintf('%1.0e',edges{3,3}),'HorizontalAlignment','center'); hold on;
text(1.25,1.3,sprintf('%2.1f',edges{4,3}),'HorizontalAlignment','center'); hold on;
text(3,0.8,sprintf('%1.0e',edges{5,3}),'HorizontalAlignment','center'); hold on;
text(3.25,1.35,sprintf('%1.0e',edges{7,3}),'HorizontalAlignment','center'); hold on;
text(4.25,1.5,sprintf('%1.0e',edges{6,3}),'HorizontalAlignment','center'); hold on;

% text(1,2.6,'5','HorizontalAlignment','center'); hold on;
% text(2,2.6,'10','HorizontalAlignment','center'); hold on;
% text(2,0.6,'0.5','HorizontalAlignment','center'); hold on;
% text(3,2.6,'2.5','HorizontalAlignment','center'); hold on;
% text(4,2.6,'1','HorizontalAlignment','center'); hold on;
% text(4,0.6,'1','HorizontalAlignment','center'); hold on;

%Text
% text(1,2.6,'NPP','HorizontalAlignment','center'); hold on;
% text(2,2.6,'Z','HorizontalAlignment','center'); hold on;
% text(2,0.6,'B','HorizontalAlignment','center'); hold on;
% text(3,2.6,'F','HorizontalAlignment','center'); hold on;
% text(4,2.6,'P','HorizontalAlignment','center'); hold on;
% text(4,0.6,'D','HorizontalAlignment','center'); hold on;


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
print(gcf, '-dpng',[fpath 'foodweb_cbrt_eg_names_v3']);
