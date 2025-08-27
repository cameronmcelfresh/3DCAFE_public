%% Script to model solidification due to a trailing beam path
gridSize = 100; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion

%% Material Parameters
Temp = 950;%Temperature in Kelvins
z = 1e14; % pre-exponential
K1 = 3e2; % pre-exponential
G_nucleation = 0.14; %Nucleation energy eV/K
k = 8.617*10^-5; % Boltzmann constant eV/K
Tm=1400; %Melting Temp, [K]
C = 7.5e-1; % pre-exponential (temperature independent)
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
scanSpeed = 0.200; %speed that the beam rasters [m/s]
laserPower = 100; %laser power [W]
crossHatch = 30e-6; %amount of crossover between scan paths [%]
layerThickness = 30e-6; %layer thickness during AM [m]
excessLength = 50e-6; %distance to start outside of the grid [m]
scanRotation=90; % degrees to rotate after each rotation [degrees]
startingZ = 80e-6; %height to start the first layer [m]
numGrains = 200; %number of grains to initialize in the voronoi tesselation

%Generate the temperature distribution given the scan speed, laser power,
%and melting temperature
%tempData = heat_fluxFunc(scanSpeed*1000,laserPower,Tm);
tempData = heat_fluxFunc_sphere(scanSpeed*1000,laserPower,Tm);
tempDist = tempData.Temps;

beamDiameter = tempData.poolDiameter; % diameter of the beam
poolRadius = round(beamDiameter/(2*realGridSize)*gridSize); %pool radius in grid size units
tDelay = 0.00001; % amount of time to wait between each beam raster %0.001
dt = 0.0002; %timestep taken each iteration , dt = 0.00005;

%beam path matrix structure:
%time[s] locationX[absolute loc] locationY[absolute loc] locationZ[absolute
%unitNormal [x] unitNormal [y]
beamPath = generateScanPath_3D(scanSpeed, crossHatch, layerThickness, scanRotation, realGridSize, dt, excessLength, startingZ);

%Movie Variables
plotData=1; %binary indicator to plot the growing microstructure or not

if plotData==1
    figure
end


%% Start Computation (all user variables set above)

%Grid to hold the positions
%each x,y position will either be 1 (indicating empty or uncrystallizated space) or a number which indicates the nuclei it is a part of.
grid = ones([gridSize,gridSize,gridSize]); 
gridSolidTime = zeros([gridSize,gridSize,gridSize]); %time that the pixel was solidified - for growth of the adjacent grains

nuclei = []; 
%each row will consist of a unique nuclei
%1: x center position
%2: y center position
%3: z center position  <--- everything was shifted over due to this...
%4: time of initilization of nuclei 
%5: color/number assignment of nuclei / grainID
%6: binary indicator of whether the nuclei is still growing or not. 1==growing 0==stopped growing, totally impinged

%Pre-compuate all possible boundary points

%Grow the original voronoi tesselation and shift it
[grid,grainCenters] = growEquiaxed3D_2(gridSize, numGrains);

grid=grid+1000; %shift it so the new grains can start counting at 1

grainCount=1;
iter=1;

dTime=toc;
timerSaved.initialize = timerSaved.initialize+dTime;

for t = 0:dt:max(beamPath)

    currentZHeight = round(beamPath(iter,4)*gridSize/realGridSize); %current z height for the depositted layer

    tic

    %Find the current temperature distribution
    currentTemps = findTempDist3D_ellipse(iter,realGridSize,gridSize,beamPath,tempData);

    %Print the temperature distribution to save
    xcolumn =xgrid(:);
    ycolumn =ygrid(:);
    zcolumn =zgrid(:);
    tempColumn = currentTemps(:);
    binaryIndicatory = zcolumn*realGridSize/gridSize<=beamPath(iter,4); %only print the temperatures at or below the layer height

    %Save temperatures in vtk format
%     vtkwrite(sprintf("temps_%i.vtk",iter),'unstructured_grid',xcolumn(binaryIndicatory),...
%         ycolumn(binaryIndicatory),...
%         zcolumn(binaryIndicatory),...
%         'scalars',...
%         'temperatures',tempColumn(binaryIndicatory));
     
    %Use the current temperature profile to melt the microstructure
    lastGrid=grid;
    if t==0
       grid(currentTemps>Tm & zgrid<=currentZHeight)=-1; %move the beam
    else
       grid(currentTemps>Tm & zgrid<=currentZHeight)=-1; %move the beam
             
       grid(currentTemps<Tm & ...
           grid<1 & zgrid<=currentZHeight)=0;
    end

    %Save microstructure in vtk format
    gridIndicator = grid>0;
    binaryGridIndicator = gridIndicator(:) & binaryIndicatory; %find the points that are within the region of interst and also solid
    
    vtkwrite(sprintf("micro_%i.vtk",iter),'unstructured_grid',xcolumn(binaryGridIndicator),...
        ycolumn(binaryGridIndicator),...
        zcolumn(binaryGridIndicator),...
        'scalars',...
        'grain_ID',grid(binaryGridIndicator));

    dTime=toc;
    timerSaved.moveBeam = timerSaved.moveBeam+dTime;

    %Grow the grains in 3D
    %grid=evolveGrid_3D(grid, currentZHeight); %growth without a perscribed growth rate
    [grid,gridSolidTime] = evolveGrid_3D_growthRate(grid,currentZHeight,gridSolidTime, realGridSize/gridSize, t, growthRate);
     
    iter=iter+1;
end

% %% Plot timers
% 
% if plotData==1
%     X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Nuclei','Grow Boundaries','Plot'});
%     Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.growBoundaries,timerSaved.plot];
%     bar(X,Y);
%     ylabel("Time [s]");
% end
% 
