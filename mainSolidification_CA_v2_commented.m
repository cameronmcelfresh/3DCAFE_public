%% Script to model solidification due to a trailing beam path
startup_mtex; %start mtex

gridSize = 100; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=200e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion
dx = realGridSize/gridSize; %grid step size
%% Material Parameters
Tm = 1400; %Melting temperature in K

%% Timing Constants
timerSaved.initialize=0;
timerSaved.moveBeam=0;
timerSaved.findNuclei=0;
timerSaved.growNuclei=0;
timerSaved.findNewOctahedrons=0;
timerSaved.removeNuclei=0;
tic

%% Scan Parameters 
scanSpeed = 0.400; %speed that the beam rasters [m/s]
laserPower = 100; %laser power [W]
crossHatch = 100e-6; %amount to step laterally each scan path [m]
layerThickness = 40e-6; %layer thickness during AM [m]
excessLength = 80e-6; %distance to start outside of the grid [m]
scanRotation=0; % degrees to rotate after each rotation [degrees]
startingZ = 80e-6; %height to start the first layer [m]
numGrains = 1500; %number of grains to initialize in the voronoi tesselation

%% Generate the temperature distribution given the scan speed, laser power, and melting temperature

[Xmesh,Ymesh,Zmesh,tempMatrix] = tempDistPreProcess("caseA_fit_temps.csv",[600,0,100],1e-6,800e-6,10e-6);
tempData.Temps = tempMatrix;
tempData.dx = 800e-6 /(length(tempMatrix));
tempData.gridSize=length(tempMatrix);
tempData.poolDiameter=140e-6;
maxTemps = zeros([gridSize,gridSize,gridSize]); %matrix to hold the maximum temperature each point sees
%tempData = heat_fluxFunc_sphere(scanSpeed*1000,laserPower,Tm); %sphereical approach
%tempDist = tempData.Temps;

%%
beamDiameter = tempData.poolDiameter; % diameter of the beam
poolRadius = round(beamDiameter/(2*realGridSize)*gridSize); %pool radius in grid size units
tDelay = 0.0002; % amount of time to wait between each beam raster %0.001
dt = 0.0000025; %timestep taken each iteration , dt = 0.00005;

%beam path matrix structure:
%time[s] locationX[absolute loc] locationY[absolute loc] locationZ[absolute
%unitNormal [x] unitNormal [y]
beamPath = generateScanPath_3D(scanSpeed, crossHatch, layerThickness, scanRotation, realGridSize, dt, excessLength, startingZ);

%% Start Computation (all user variables set above)

nuclei = []; 
%each row will consist of a unique nuclei
% column 1-3: x,y,z position of the nuclei
% column 4: index of the cell that the nuclei belongs to
% column 5-7: euler angles of the nuclei
% column 8: initial nucleation time of the nuclei
% column 9: size of the nuclei
% column 10: grainID of nuclei

Cinc = 2e-3; %incubation time constant, K*s
incubationTheta = zeros(size(xgrid)); %matrix to hold the incubation time of each liquid point for mushy zone nucation

%Grow the original voronoi tesselation and shift it
%[grid,grainCenters] = growEquiaxed3D_2(gridSize, numGrains);
%grid = growEquiaxed3D_4(gridSize, numGrains,round(layerThickness/realGridSize*gridSize),round(excessLength/realGridSize*gridSize)); %using finite layers
grid = growEquiaxed3D_5(gridSize, numGrains,round(layerThickness/realGridSize*gridSize),round(excessLength/realGridSize*gridSize)); %using finite layers

grid=grid+1000; %shift it so the new grains can start counting at 1

%Build a list of Euler angles for each unique grain
uniqueGrains = unique(grid);
rot = rotation.rand(length(uniqueGrains));
%euler = [uniqueGrains,rand([length(uniqueGrains),3])]; %construct a collection of random orientations
euler = [uniqueGrains,rot.phi1,rot.Phi,rot.phi2]; %construct a collection of random orientations using MTEX
rotMatrices = rot.matrix'; %rotation matrices to handle the octahedron rotations
iter=1;

dTime=toc;
timerSaved.initialize = timerSaved.initialize+dTime;

for t = 0:dt:max(beamPath)

    %Find the current z height for the depositted layer
    tic
    currentZHeight = round(beamPath(iter,4)*gridSize/realGridSize); 

    %% Find the current temperature distribution
    %currentTemps = findTempDist3D_ellipse(iter,realGridSize,gridSize,beamPath,tempData); %sphereical approach
    %currentTemps = findTempDist3D(iter,realGridSize,gridSize,beamPath,tempData); %FEA approach
    currentTemps = findTempDist3D_v1(iter,currentZHeight,realGridSize,gridSize,beamPath,tempData); %FEA approach
    maxTemps = max(currentTemps,maxTemps); %update the maximum temperature matrix
    
    %% Save the temperature distribution in vtk format 
%     xcolumn =xgrid(:);
%     ycolumn =ygrid(:);
%     zcolumn =zgrid(:);
%     tempColumn = currentTemps(:);
%     binaryIndicatory = zcolumn*realGridSize/gridSize<=beamPath(iter,4); %only print the temperatures at or below the layer height
% 
%      vtkwrite(sprintf("temps_%i.vtk",iter),'unstructured_grid',xcolumn(binaryIndicatory),...
%          ycolumn(binaryIndicatory),...
%          zcolumn(binaryIndicatory),...
%          'scalars',...
%          'temperatures',tempColumn(binaryIndicatory));
     
    %% Use the current temperature profile to melt the microstructure
    lastGrid=grid;
    if t==0
       grid(currentTemps>Tm & zgrid<=currentZHeight)=-1; %move the beam
    else
       grid(currentTemps>Tm & zgrid<=currentZHeight)=-1; %move the beam
             
       grid(currentTemps<Tm & ...
           grid<1 & zgrid<=currentZHeight)=0;
    end

    dTime = toc;
    timerSaved.moveBeam=timerSaved.moveBeam+dTime;


    %% Remove any surrounded nuclei
    tic
    nuclei = removeSurroundedNuclei_v1(nuclei,grid);
    dTime = toc;
    timerSaved.removeNuclei=timerSaved.removeNuclei+dTime;

    %% Find new octohedrons from growth
    tic
    [grid,nuclei] = findNewOctohedrons_v1(grid,nuclei,t,currentZHeight);
    dTime = toc;
    timerSaved.findNewOctahedrons=timerSaved.findNewOctahedrons+dTime;

    %% Nucleate new octahedrons
    tic
    incubationTheta = incubationTheta + dt./(Cinc./(Tm-currentTemps)); %update the incubation time matrix
    incubationTheta(incubationTheta>0)=0; %reset incubation theta of solid cells to zero
    incubationTheta(currentTemps>Tm)=0; %reset incubation theta of liquid cells to zero
    [grid,nuclei,euler] = nucleateOctahedrons_v3(grid,nuclei,euler,t,currentZHeight,currentTemps,dx,dt,maxTemps,incubationTheta);
    dTime = toc;
    timerSaved.findNuclei=timerSaved.findNuclei+dTime;

    %% Grow the grains in 3D
    tic
    [nuclei] = growNuclei_v1(nuclei,currentTemps,Tm,dt,gridSize,realGridSize); %local temperature polynomial growth velocity approach
    dTime = toc;
    timerSaved.growNuclei=timerSaved.growNuclei+dTime;
    %% Save microstructure in vtk format
%     gridIndicator = grid>0;
%     binaryGridIndicator = gridIndicator(:) & binaryIndicatory; %find the points that are within the region of interst and also solid
%     
%     vtkwrite(sprintf("micro_%i.vtk",iter),'unstructured_grid',xcolumn(binaryGridIndicator),...
%         ycolumn(binaryGridIndicator),...
%         zcolumn(binaryGridIndicator),...
%         'scalars',...
%         'grain_ID',grid(binaryGridIndicator));

    fprintf("Iteration %i\n",iter);

    iter=iter+1;
end


%% Print the final structures
iter=iter-1;
xcolumn =xgrid(:);
ycolumn =ygrid(:);
zcolumn =zgrid(:);
tempColumn = currentTemps(:);
binaryIndicatory = zcolumn*realGridSize/gridSize<=beamPath(iter,4); %only print the temperatures at or below the layer height

 vtkwrite(sprintf("temps_%i.vtk",iter),'unstructured_grid',xcolumn(binaryIndicatory),...
     ycolumn(binaryIndicatory),...
     zcolumn(binaryIndicatory),...
     'scalars',...
     'temperatures',tempColumn(binaryIndicatory));

% To randomly permute all the grain IDs
gridIndicator = grid>0;
binaryGridIndicator = gridIndicator(:) & binaryIndicatory; %find the points that are within the region of interst and also solid
gridCopy = grid;
uniqueID = unique(grid);
randID  = randperm(length(uniqueID));
for i = 1:length(uniqueID)
    gridCopy(gridCopy==uniqueID(i)) = randID(i);
end

vtkwrite(sprintf("micro_%i.vtk",iter),'unstructured_grid',xcolumn(binaryGridIndicator),...
    ycolumn(binaryGridIndicator),...
    zcolumn(binaryGridIndicator),...
    'scalars',...
    'grain_ID',gridCopy(binaryGridIndicator));


%% Plot timers
figure
X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Octahedrons','Find New Octahedrons','Remove Octahedrons'});
Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.findNewOctahedrons,timerSaved.removeNuclei];
bar(X,Y);
ylabel("Time [s]");

%% Plot Pole Figures

cs = crystalSymmetry('cubic');
%ori = orientaion.rand(cs);
rot = rotation.byEuler(euler(:,2),euler(:,3),euler(:,4));
h = Miller({1,0,0},cs);
r = rot * h.symmetrise;
figure
plot(r)
title("All grains")

g = unique(grid(:));

%rot = rotation.byEuler(180*3.1415*gridEuler1,180*3.1415*gridEuler2,180*3.1415*gridEuler3);
rot = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
r = rot * h.symmetrise;
figure
plot(r)
title("100 - PF, Figure,Original",'contourf')
%mtexColorbar

gridEuler1 = grid(xgrid>40 & xgrid<75);
g = unique(gridEuler1(:));

%rot_mid = rotation.byEuler(180*3.1415*gridEuler1,180*3.1415*gridEuler2,180*3.1415*gridEuler3);
rot_mid = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);

r_mid = rot_mid * h.symmetrise;
figure
plot(r_mid)
title("100 - PF,Mid Slice",'contourf')


%% Inverse Pole Figures

g = unique(grid(:));

mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);
figure
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF, Original")

figure
plotIPDF(odf,vector3d(0,1,1),'antipodal')
title("011 - IPF, Original")


%% Middle slice, should have texture

gridEuler1 = grid(xgrid>40 & xgrid<75);
g = unique(gridEuler1(:));

mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);

figure
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF, Middle")

figure
plotIPDF(odf,vector3d(0,1,1),'antipodal')
title("011 - IPF, Middle")


%% Bottom slice - should be untextured?

gridEuler1 = grid(xgrid<15);
g = unique(gridEuler1(:));

mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);

figure
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF, Bottom")

figure
plotIPDF(odf,vector3d(0,1,1),'antipodal')
title("011 - IPF, Bottom")


%% Random test

randMod1 = rotation.rand(1000);
mod1 = rotation.byEuler(randMod1.phi1,randMod1.Phi,randMod1.phi2,cs);
odf = unimodalODF(mod1);
figure
plotIPDF(odf,[zvector],'antipodal')
title("100 - IPF, Middle")
