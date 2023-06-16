function [] = deformedDensityPlot(DT,D,eMat)

NC = DT.Points;                     %Nodal Coordinates (NC); n_point x 3
LI = edges(DT);                     %List of Edges

% Plot Results
Contours = 24;
Scolor = floor(rescale([0; eMat; 1;],1,Contours));
Scolor([1 end]) = [];
C = brewermap(Contours,'Greys');
colormap(C);

x = NC(:,1)+D(:,1);
y = NC(:,2)+D(:,2);

for i = 1:length(LI)
   %plot deformed elements
   plot([x(LI(i,1)) x(LI(i,2))],[y(LI(i,1)) y(LI(i,2))],...
       'linewidth',2,'color',C(Scolor(i),:)); hold on
end
axis equal

end