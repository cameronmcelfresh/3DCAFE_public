%% Script to model solidification due to a trailing beam path

foldername= 'test_2';

gridSize = 200; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid. Relevant to the velocity of the boundary motion

% Start Computation (all user variables set above)

%Grid to hold the positions
%each x,y position will either be 1 (indicating empty or uncrystallizated space) or a number which indicates the nuclei it is a part of.
grid = ones([gridSize,gridSize]); 

%Plotting info
cmap = hsv(2000); 
cmap = cmap(randperm(2000),:);

cmap(1,:)=[1 1 1]; %set the laser color to white
cmap(2,:)=[0 0 0]; %set the liquid color to black
cmap(3,:)=[1 1 1]; %set the outside color to white

cmapTemp = parula(300);

%Grow the original voronoi tesselation and shift it
[grid,grainCenters] = growEquiaxed(gridSize,100);
grid=grid+5; %shift it so the new grains can start counting at 1

writeDAMASK(2,grid,xgrid,ygrid,200,realGridSize/gridSize,cmap,foldername)