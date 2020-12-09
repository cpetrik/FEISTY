
% three scenearios based on the proportion of the zooplankton community that migrates to depth 
% see calcpreference (3x)
%% scneario 1 (calcpreference_visual_1)
 % zp_d = (zp_n + zp_n + zp_d)/3;
 %% scneario 2 (calcpreference -- main document)
 % zp_d = (zp_n + zp_d)/2;
 %% scneario 3 (calcpreference_visual_3)
 % zp_d = (zp_n + zp_d + zp_d)/3;

% data input
bRp = linspace(5,50,6); % change maximum pelagic productivity
visual = linspace(1,1.9,10); % change how good/bad visual predators see in the dark
 
 %% scenario 1
rowle = length(visual);
colle = length(bRp); 
BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);

for j = 1:colle
for i = 1:rowle
    baseparameters % region is 4, photic depth is 150
    param.bottom = 1000;  % depth in meter

    % calculate resource productivities
    param.K =  [bRp(j)  bRp(j)   0    0];  % g C ww/m2/yr
    
    % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude all other fish
    param.y0(param.ix1(3):param.ix2(3))=0; %large pelagics to zero
    param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    param.y0(param.ix1(5):param.ix2(5))=0; %demersal to zero
    
    % feeding preference matrix
    param.visual = visual(i);
    [param.theta, param.depthDay, param.depthNight] = calcpreference_visual_1(param); 
   
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i,j) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2))));
end
end

save visual_one


 %% scenario 2
rowle = length(visual);
colle = length(bRp); 
BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);

for j = 1:colle
for i = 1:rowle
    baseparameters % region is 4, photic depth is 400
    param.bottom = 1000;  % depth in meter

    % calculate resource productivities
    param.K =  [bRp(j)  bRp(j)   0    0];  % g C ww/m2/yr
    
    % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude all other fish
    param.y0(param.ix1(3):param.ix2(3))=0; %large pelagics to zero
    param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    param.y0(param.ix1(5):param.ix2(5))=0; %demersal to zero
    
    % feeding preference matrix
    param.visual = visual(i);
    [param.theta, param.depthDay, param.depthNight] = calcpreference(param); 
   
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i,j) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2)))); 
end
end

save visual_two

%% scenario 3
rowle = length(visual);
colle = length(bRp); 
BioFf = zeros(rowle,colle);
BioMf = zeros(rowle,colle);

for j = 1:colle
for i = 1:rowle
    baseparameters % region is 4, photic depth is 400
    param.bottom = 1000;  % depth in meter

    % calculate resource productivities
    param.K =  [bRp(j)  bRp(j)   0    0];  % g C ww/m2/yr
    
    % starting values
    if (i==1)
        % Default initial conditions:
        param.y0 = [0.1*param.K param.B0];
    else
        % Use "results" argument for initial conditions:
        param.y0 = result.y(end,:) + 0.01;
    end
    
    % exclude all other fish
    param.y0(param.ix1(3):param.ix2(3))=0; %large pelagics to zero
    param.y0(param.ix1(4):param.ix2(4))=0; %bathypelagics to zero
    param.y0(param.ix1(5):param.ix2(5))=0; %demersal to zero
    
    % feeding preference matrix
    param.visual = visual(i);
    [param.theta, param.depthDay, param.depthNight] = calcpreference_visual_3(param); 
   
    result = poem(param);
    y = result.y;
    Bin = floor(0.8*length(y));
       
    BioFf(i,j) = sum(mean(y(Bin:end,param.ix1(1):param.ix2(1))));
    BioMf(i,j) = sum(mean(y(Bin:end,param.ix1(2):param.ix2(2))));
end
end

save visual_three

%% plot
load('visual_one.mat')

FigH = figure;
SP1H = subplot(3,2,1);
imagesc(BioFf)
hold on
yticks([1 5 10])
yticklabels([visual(1) visual(5) visual(10)])
xticks(linspace(1,6,6))
xticklabels(bRp)
colorbar('off')
caxis([0 8])
hold off
c = colorbar;
c.Label.String = 'Biomass (g WW m-2)';

SP2H = subplot(3,2,2);
imagesc(BioMf)
hold on
yticks(linspace(1,10,10))
yticks([1 5 10])
yticklabels([visual(1) visual(5) visual(10)])
xticklabels(bRp)
caxis([0 8])
colorbar('off')
SP2H.Position(3)= SP1H.Position(3);
hold off

load('visual_two.mat')
SP3H = subplot(3,2,3);
imagesc(BioFf)
hold on
yticks([1 5 10])
yticklabels([visual(1) visual(5) visual(10)])
xticks(linspace(1,6,6))
xticklabels(bRp)
colorbar('off')
caxis([0 8])
SP3H.Position(3)= SP1H.Position(3);
hold off

SP4H = subplot(3,2,4);
imagesc(BioMf)
hold on
yticks([1 5 10])
yticklabels([visual(1) visual(5) visual(10)])
xticks(linspace(1,6,6))
xticklabels(bRp)
caxis([0 8])
colorbar('off')
SP4H.Position(3)= SP1H.Position(3);
hold off

load('visual_three.mat')
SP5H = subplot(3,2,5);
imagesc(BioFf)
hold on
yticks([1 5 10])
yticklabels([visual(1) visual(5) visual(10)])
xticks(linspace(1,6,6))
xticklabels(bRp)
colorbar('off')
caxis([0 8])
SP5H.Position(3)= SP1H.Position(3);
hold off

SP6H = subplot(3,2,6);
imagesc(BioMf)
hold on
yticks([1 5 10])
yticklabels([visual(1) visual(5) visual(10)])
xticks(linspace(1,6,6))
xticklabels(bRp)
caxis([0 8])
colorbar('off')
SP6H.Position(3)= SP1H.Position(3);
hold off

YLabel3H = get(SP3H,'YLabel');
set(YLabel3H,'String','Scalar of visual predation (sc)');
set(YLabel3H,'Position',[-.4 5 0]);
XLabel5H = get(SP5H,'XLabel');
set(XLabel5H,'String','Max. zooplankton prod. (g WW m-2 y-1)');
set(XLabel5H,'Position',[7 13 0]);


