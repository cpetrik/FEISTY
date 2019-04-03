N = 1000 ;
x = 2*rand(N,1) - 1; 
y = 2*rand(N,1) - 1; 
z = 2*rand(N,1) - 1;
v = x.^2 + y.^3 - z.^4;
%
d = -1:0.05:1;
[xq,yq,zq] = meshgrid(d,d,0);
% Interpolate the scattered data on the grid. Plot the results.
vq = griddata(x,y,z,v,xq,yq,zq);
plot3(x,y,v,'ro')
hold on
surf(xq,yq,vq)

%% doesn't work
%test = fx(idbot,:);
%test(:,6) = baic_all(idbot);
x = test(:,1); 
y = test(:,3); 
z = test(:,5); 
v = test(:,6); 

lam = [0.6,0.75];
be = [0.15,0.25];
ae = [50,100];
%[agrid,bgrid,egrid] = meshgrid([0.6,0.75,0.9], [0.15,0.25,0.35], [50,100,150]);
%cgrid = nan*ones(size(agrid));

[xq,yq,zq] = meshgrid(lam,be,ae);
% Interpolate the scattered data on the grid. Plot the results.
vq = griddata(x,y,z,v,xq,yq,zq);
plot3(x,y,v,'ro')
hold on
surf(xq,yq,zq,vq)

%%
[y,z] = meshgrid(linspace(0,10,40));
for off=50:50:200    
  x = off + zeros(size(z));
    % My standin for your vorticity data
    c = cos((x+y)/5) .* cos((x+z)/5);
    surf(x,y,z,c)
    hold on
end
hold off
xlim([0 200])

%%
figure
surf(squeeze(xq(:,1,:)),squeeze(yq(:,1,:)),squeeze(zq(:,1,:)),squeeze(vq(:,1,:)));
hold on;
surf(squeeze(xq(:,2,:)),squeeze(yq(:,2,:)),squeeze(zq(:,2,:)),squeeze(vq(:,2,:)));

%%
agrid = NaN(3,3,3);
bgrid = NaN(3,3,3);
cgrid = NaN(3,3,3);
egrid = NaN(3,3,3);

agrid(1:2,1:2,1:2) = xq;
bgrid(1:2,1:2,1:2) = yq;
egrid(1:2,1:2,1:2) = zq;
cgrid(1:2,1:2,1:2) = vq;

%%
figure
surf(squeeze(agrid(:,1,:)),squeeze(bgrid(:,1,:)),squeeze(egrid(:,1,:)),...
    squeeze(cgrid(:,1,:)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(:,2,:)),squeeze(bgrid(:,2,:)),squeeze(egrid(:,2,:)),...
    squeeze(cgrid(:,2,:)),'FaceColor','interp');
colormap('jet')
caxis([360 400])

figure
surf(squeeze(agrid(1,:,:)),squeeze(bgrid(1,:,:)),squeeze(egrid(1,:,:)),...
    squeeze(cgrid(1,:,:)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(2,:,:)),squeeze(bgrid(2,:,:)),squeeze(egrid(2,:,:)),...
    squeeze(cgrid(2,:,:)),'FaceColor','interp');
colormap('jet')
caxis([360 400])

figure
surf(squeeze(agrid(:,:,1)),squeeze(bgrid(:,:,1)),squeeze(egrid(:,:,1)),...
    squeeze(cgrid(:,:,1)),'FaceColor','interp');
hold on;
surf(squeeze(agrid(:,:,2)),squeeze(bgrid(:,:,2)),squeeze(egrid(:,:,2)),...
    squeeze(cgrid(:,:,2)),'FaceColor','interp');
colormap('jet')
caxis([360 400])

%%
figure
surf(squeeze(agrid(:,1,:)),squeeze(bgrid(:,1,:)),squeeze(egrid(:,1,:)),...
    squeeze(cgrid(:,1,:)),'FaceColor','interp');
hold on;

surf(squeeze(agrid(1,:,:)),squeeze(bgrid(1,:,:)),squeeze(egrid(1,:,:)),...
    squeeze(cgrid(1,:,:)),'FaceColor','interp');
hold on;
% surf(squeeze(agrid(2,:,:)),squeeze(bgrid(2,:,:)),squeeze(egrid(2,:,:)),...
%     squeeze(cgrid(2,:,:)),'FaceColor','interp');

% surf(squeeze(agrid(:,:,1)),squeeze(bgrid(:,:,1)),squeeze(egrid(:,:,1)),...
%     squeeze(cgrid(:,:,1)),'FaceColor','interp');
% hold on;
surf(squeeze(agrid(:,:,2)),squeeze(bgrid(:,:,2)),squeeze(egrid(:,:,2)),...
    squeeze(cgrid(:,:,2)),'FaceColor','interp');
colormap('jet')
caxis([360 400])

%%
figure
% surf(squeeze(agrid(:,1,:)),squeeze(bgrid(:,1,:)),squeeze(egrid(:,1,:)),...
%     squeeze(cgrid(:,1,:)),'FaceColor','interp');
surf(squeeze(agrid(:,2,:)),squeeze(bgrid(:,2,:)),squeeze(egrid(:,2,:)),...
    squeeze(cgrid(:,2,:)),'FaceColor','interp');
hold on;

% surf(squeeze(agrid(1,:,:)),squeeze(bgrid(1,:,:)),squeeze(egrid(1,:,:)),...
%     squeeze(cgrid(1,:,:)),'FaceColor','interp');
% hold on;
surf(squeeze(agrid(2,:,:)),squeeze(bgrid(2,:,:)),squeeze(egrid(2,:,:)),...
    squeeze(cgrid(2,:,:)),'FaceColor','interp');
hold on;

surf(squeeze(agrid(:,:,1)),squeeze(bgrid(:,:,1)),squeeze(egrid(:,:,1)),...
    squeeze(cgrid(:,:,1)),'FaceColor','interp');
hold on;
% surf(squeeze(agrid(:,:,2)),squeeze(bgrid(:,:,2)),squeeze(egrid(:,:,2)),...
%     squeeze(cgrid(:,:,2)),'FaceColor','interp');
colormap('jet')
caxis([360 400])



