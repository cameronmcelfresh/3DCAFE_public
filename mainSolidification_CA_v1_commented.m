%% Script to model solidification due to a trailing beam path
gridSize = 100; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion

%% Material Parameters
Temp = 1050;%Temperature in Kelvins
z = 1e14; % pre-exponential
K1 = 3e2; % pre-exponential
G_nucleation = 0.14; %Nucleation energy eV/K
k = 8.617*10^-5; % Boltzmann constant eV/K
Tm=1400; %Melting Temp, [K]
C = 7.5e6; % pre-exponential (temperature independent)
Qgrowth = 0.6; % activation energy for boundary motion/growth eV/K

Gstar = K1*1/(Tm-Temp)^2; %driving force for nucleation due to undercooling
Nhomo=z*exp(-Gstar/(Temp*k))*exp(-G_nucleation/(Temp*k)); %homogeneous nucleation rate
growthRate = C*exp(-Qgrowth/(k*Temp));

nucleationRate = Nhomo; %per unit time
nucleationRateSTD = 0.1; %nucleation rate STD to use during sampling
shapeFactor = 0.2; %heterogeneous nucleation shape factor. If set to anything above 1 then heterogeneous nucleation is suspended

%Save the nucleation info
nucleationData.z = z;
nucleationData.K1 = K1;
nucleationData.G_nucleation = G_nucleation;
nucleationData.k = k;
nucleationData.Tm = Tm;
nucleationData.nucRateSTD = nucleationRateSTD;
nucleationData.shapeFactor=shapeFactor;

% Timing Constants
timerSaved.initialize=0;
timerSaved.moveBeam=0;
timerSaved.findNuclei=0;
timerSaved.growNuclei=0;
timerSaved.growBoundaries=0;
timerSaved.plot =0;
tic

%% Scan Parameters 
scanSpeed = 0.400; %speed that the beam rasters [m/s]
laserPower = 100; %laser power [W]
crossHatch = 60e-6; %amount to step laterally each scan path [m]
layerThickness = 40e-6; %layer thickness during AM [m]
excessLength = 180e-6; %distance to start outside of the grid [m]
scanRotation=67; % degrees to rotate after each rotation [degrees]
startingZ = 80e-6; %height to start the first layer [m]
numGrains = 1600; %number of grains to initialize in the voronoi tesselation

%% Generate the temperature distribution given the scan speed, laser power, and melting temperature

[Xmesh,Ymesh,Zmesh,tempMatrix] = tempDistPreProcess("caseA_fit_temps.csv",[600,0,100],1e-6,800e-6,10e-6);
tempData.Temps = tempMatrix;
tempData.dx = 800e-6 /(length(tempMatrix));
tempData.gridSize=length(tempMatrix);
tempData.poolDiameter=140e-6;

%tempData = heat_fluxFunc_sphere(scanSpeed*1000,laserPower,Tm); %sphereical approach
%tempDist = tempData.Temps;

%%

beamDiameter = tempData.poolDiameter; % diameter of the beam
poolRadius = round(beamDiameter/(2*realGridSize)*gridSize); %pool radius in grid size units
tDelay = 0.0002; % amount of time to wait between each beam raster %0.001
dt = 0.000001; %timestep taken each iteration , dt = 0.00005;

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

%Grow the original voronoi tesselation and shift it
%[grid,grainCenters] = growEquiaxed3D_2(gridSize, numGrains);
%grid = growEquiaxed3D_3(gridSize, numGrains,round(realGridSize/layerThickness)); %using finite layers
grid = growEquiaxed3D_4(gridSize, numGrains,round(layerThickness/realGridSize*gridSize),round(excessLength/realGridSize*gridSize)); %using finite layers

grid=grid+1000; %shift it so the new grains can start counting at 1

%Build a list of Euler angles for each unique grain
uniqueGrains = unique(grid);
euler = [uniqueGrains,rand([length(uniqueGrains),3])]; %construct a collection of random orientations

iter=1;

dTime=toc;
timerSaved.initialize = timerSaved.initialize+dTime;

for t = 0:dt:max(beamPath)

    %Find the current z height for the depositted layer
    currentZHeight = round(beamPath(iter,4)*gridSize/realGridSize); 

    %% Find the current temperature distribution
    %currentTemps = findTempDist3D_ellipse(iter,realGridSize,gridSize,beamPath,tempData); %sphereical approach
    %currentTemps = findTempDist3D(iter,realGridSize,gridSize,beamPath,tempData); %FEA approach
    currentTemps = findTempDist3D_v1(iter,currentZHeight,realGridSize,gridSize,beamPath,tempData); %FEA approach

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

    %% Remove any surrounded nuclei
    nuclei = removeSurroundedNuclei_v1(nuclei,grid);

    %% Find new grains/octohedrons from growth
    [grid,nuclei] = findNewOctohedrons_v1(grid,nuclei,t,currentZHeight);

    %% Nucleate new octahedrons
    %[grid,nuclei] = nucleateOctahedrons(grid,nuclei,euler,t,currentZHeight);
    [grid,nuclei] = nucleateOctahedrons_v1(grid,nuclei,euler,t,currentZHeight);

    %% Grow the grains in 3D
    %[nuclei] = growNuclei(nuclei,growthRate,dt); %constant growth velocity
    %approach
    [nuclei] = growNuclei_v1(nuclei,currentTemps,Tm,dt,gridSize,realGridSize); %local temperature polynomial growth velocity approach

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

gridIndicator = grid>0;
binaryGridIndicator = gridIndicator(:) & binaryIndicatory; %find the points that are within the region of interst and also solid

vtkwrite(sprintf("micro_%i.vtk",iter),'unstructured_grid',xcolumn(binaryGridIndicator),...
    ycolumn(binaryGridIndicator),...
    zcolumn(binaryGridIndicator),...
    'scalars',...
    'grain_ID',grid(binaryGridIndicator));

toc
% %% Plot timers
% 
% if plotData==1
%     X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Nuclei','Grow Boundaries','Plot'});
%     Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.growBoundaries,timerSaved.plot];
%     bar(X,Y);
%     ylabel("Time [s]");
% end
% 
