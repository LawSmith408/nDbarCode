function [] = rhoPlot(DT,D,rho)

NC = DT.Points;                     %Nodal Coordinates (NC); n_point x 3
LI = edges(DT);                     %List of Edges

% Plot Results
% Contours = 24;
% Scolor = floor([0; eMat; 1;]);
% Scolor = floor(rescale([0; eMat; 1;],1,Contours));
% Scolor([1 end]) = [];
C = brewermap(Contours,'Greys');
colormap(C);

% %add the deformations to the nodal coordinates
% NC = NC+D;

% %add a z coordinate to the nodal coordinates
% NC = [NC zeros(size(NC,1),1)];

% %add a set of new points that are on the midpoints of each edge
% xM = mean([NC(LI(:,1),1)'; NC(LI(:,2),1)'])';
% yM = mean([NC(LI(:,1),2)'; NC(LI(:,2),2)'])';
% zM = zeros(size(xM))-0.3;

% %create a set of faces
% F = [LI (1:length(LI))'+length(NC)];

% %Augment the node list with a duplicate set of nodes;
% NC2 = [NC; [xM yM zM]];

%let's make a new node list that includes duplicates. each line segment
%will have its own two nodes

V = NC(LI(:,1),:)

%plot using patch
tic;
cFigure; axisGeom
patch('faces',F,'vertices',NC2,'edgecolor',rho); hold on
camlight headlight
fprintf('Elapsed time: %.2f seconds\n', toc);

view(3)

% tic; figure;
% % Plot the edges on top of the triangles
% for i = 1:length(LI)
%     plot([NC(LI(i,1),1) NC(LI(i,2),1)], [NC(LI(i,1),2) NC(LI(i,2),2)], ...
%         'linewidth', 2, 'color', C(Scolor(i),:));
% end


axis equal

end
