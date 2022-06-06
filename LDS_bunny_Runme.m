clear; clc; close all

%% load geometry
%load small_bunny_mesh2D.mat;
%load small_bunny_mesh2D_vfine.mat
load small_bunny_mesh3D.mat;
%load small_bunny_mesh3D_fine.mat;


%define boundary conditions
NC = DT.Points;
fixed = find(NC(:,1)<-40);
forced = find(NC(:,1)>15);

%solve system
[D,C] = LDS_Solver(DT,5e4,fixed,forced);

%% Generate Output Plot
figure
patch('faces',DT.ConnectivityList,'vertices',NC+D,'FaceColor','w','handlevisibility','off');hold on
plotV(NC(fixed,:)+D(fixed,:),'ro')
plotV(NC(forced,:)+D(forced,:),'bo')
legend('fixed','loaded'); axis equal
if size(NC,2)>2; axisGeom; end; gdrawnow()