%% Script to model solidification due to a trailing beam path
addpath '/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MTEX/mtex-5.8.1'
startup_mtex; %start mtex

params = {
%"caseA_fit_temps", 0.40,25e-6,50e-6,10^-2,2e-3;...
     %"caseB_fit_temps",0.8,20e-6,50e-6,10^-4,2e-3;
     %"caseC_fit_temps", 1.2,20e-6, 50e-6,10^-4,2e-3;
     %"caseC_fit_temps", 1.3,20e-6, 50e-6,10^-4,2e-3
%     "caseB_fit_temps",0.8,25e-6,80e-6,10,2e-2;
%     "caseC_fit_temps", 1.2,25e-6, 80e-6,10,2e-2;
%     "caseC_fit_temps", 1.3,25e-6, 90e-6,10,2e-2
    %"caseA_fit_temps", 0.40,25e-6,80e-6,10 ;...
    %"caseB_fit_temps",0.8,25e-6,80e-6,10 ;
    %"caseC_fit_temps", 1.2,25e-6, 80e-6, 10;
    %"caseC_fit_temps", 1.3,25e-6, 90e-6,10
%"caseN05_temps", 2.20,30e-6,60e-6,5*10^3,4e-2;...
%"caseN05_temps", 2.20,30e-6,85e-6,5*10^3,4e-2;...
"caseN05_temps", 2.20,30e-6,150e-6,5*10^3,4e-2;...
%"caseN05_temps", 2.20,30e-6,150e-6,5*10^3,4e-2;...
%"caseN05_temps", 2.20,25e-6,75e-6,5*10^3,4e-2;...
%"caseN04_temps", 1.470,25e-6,50e-6,5*10^3,4e-2;...
%"caseN04_temps", 1.470,30e-6,60e-6,5*10^3,4e-2;...
%"caseN04_temps", 1.470,30e-6,80e-6,5*10^3,4e-2;...
    };

for paramIter=1:length(params)

save_dir = sprintf("test_%i",paramIter);
finalGridSize = [50,50]; %Grid size to simulate with CPFE

gridSize = 50; %100c% side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion
dx = realGridSize/gridSize; %grid step size
%% Material Parameters
%Tm = 1623; %Melting temperature in K - IN625
Tm = 1658; %Melting temperature in K - 316L

%% Timing Constants
timerSaved.initialize=0;
timerSaved.moveBeam=0;
timerSaved.findNuclei=0;
timerSaved.growNuclei=0;
timerSaved.findNewOctahedrons=0;
timerSaved.removeNuclei=0;
tic

%% Scan Parameters 
scanSpeed = params{paramIter,2};  %0.400; %speed that the beam rasters [m/s]
laserPower = 100; %laser power [W]
crossHatch =  params{paramIter,4};  %80e-6; %amount to step laterally each scan path [m]
layerThickness = params{paramIter,3};% 40e-6; %layer thickness during AM [m]
excessLength = 110e-6; %distance to start outside of the grid [m], 110e-6 for IN625
scanRotation=0; % degrees to rotate after each rotation [degrees]
startingZ = 30e-6; %75e-6; %height to start the first layer [m]
numGrains = 2500; %number of grains to initialize in the voronoi tesselation

%% Generate the temperature distribution given the scan speed, laser power, and melting temperature

[Xmesh,Ymesh,Zmesh,tempMatrix] = tempDistPreProcess(params{paramIter,1},[600,0,100],1e-6,500e-6,5e-6);
tempData.dx = 500e-6 /(length(tempMatrix));

tempData.Temps = tempMatrix;

tempData.gridSize=length(tempMatrix);
tempData.poolDiameter=140e-6;
maxTemps = zeros([gridSize,gridSize,gridSize]); %matrix to hold the maximum temperature each point sees

%%
tDelay = 0.000001; % amount of time to wait between each beam raster %0.000025
dt = 0.00001; %timestep taken each iteration , dt = 0.00005;

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
% column 11-20: rotation matrix

%% Nucleation parameters
%nuc_params.Tm = 1623; % %melting temperature,K - In625
nuc_params.Tm = Tm; %Melting temperature in K - 316L

nuc_params.Nmax = 1e18*params{paramIter,5}; %maximum nucleation rate, nuclei / m^3 * s
nuc_params.underCool_N = 5; %mean nucleation undercooing
nuc_params.underCool_sigma = 1; %standard deviation nucleation undercooling

nuc_params.Nmax_b = 5e20*params{paramIter,5}; %maximum nucleation rate for the mushy zone, nuclei / m^3 * s
nuc_params.sigma_b = 50; %standard deviation for mushy zone nucleation, K

Cinc = params{paramIter,6};%2e-3; %incubation time constant, K*s
incubationTheta = zeros(size(xgrid)); %matrix to hold the incubation time of each liquid point for mushy zone nucation

%Grow the original voronoi tesselation and shift it
%grid = growEquiaxed3D_5(gridSize, numGrains,round(layerThickness/realGridSize*gridSize),round(excessLength/realGridSize*gridSize)); 
grid = growEquiaxed3D_6(gridSize, numGrains);

grid=grid+1000; %shift it so the new grains can start counting at 1

%Build a list of Euler angles for each unique grain
uniqueGrains = unique(grid);
rot = rotation.rand(length(uniqueGrains));
euler = [uniqueGrains,rot.phi1,rot.Phi,rot.phi2]; %construct a collection of random orientations using MTEX
rotMat = rot.matrix; %rotation matrices to handle the octahedron rotations. Due to MTEX Bunge rotations, the transpose of the rotation matrix is needed for proper euler configuration.
iter=1;

dTime=toc;
timerSaved.initialize = timerSaved.initialize+dTime;

%% Run simulation

for loopi = 1:length(beamPath)

    t = beamPath(iter,1)-dt;

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
%      binaryIndicatory = zcolumn*realGridSize/gridSize<=beamPath(iter,4); %only print the temperatures at or below the layer height
% 
%      vtkwrite(sprintf("temps_%i.vtk",iter),'unstructured_grid',xcolumn(binaryIndicatory),...
%          ycolumn(binaryIndicatory),...
%          zcolumn(binaryIndicatory),...
%          'scalars',...
%          'temperatures',tempColumn(binaryIndicatory));
     
    %% Use the current temperature profile to melt the microstructure
%     lastGrid=grid;
%     if t==0
%        grid(currentTemps>Tm & zgrid<=currentZHeight)=-1; %move the beam
%     else
%        grid(currentTemps>Tm & zgrid<=currentZHeight)=-1; %move the beam
%              
%        grid(currentTemps<Tm & ...
%            grid<1 & zgrid<=currentZHeight)=0;
%     end
% 
%     dTime = toc;
%     timerSaved.moveBeam=timerSaved.moveBeam+dTime;

% 
%     %% Remove any surrounded nuclei
%     tic
%     nuclei = removeSurroundedNuclei_v1(nuclei,grid);
%     dTime = toc;
%     timerSaved.removeNuclei=timerSaved.removeNuclei+dTime;
% 
%     %% Find new octohedrons from growth
%     tic
%     [grid,nuclei] = findNewOctohedrons_v2(grid,nuclei,t,currentZHeight);
%     dTime = toc;
%     timerSaved.findNewOctahedrons=timerSaved.findNewOctahedrons+dTime;
% 
%     %% Nucleate new octahedrons
%     tic
%     incubationTheta = incubationTheta + dt./(Cinc./(Tm-currentTemps)); %update the incubation time matrix
%     incubationTheta(grid>0)=0; %reset incubation theta of solid cells to zero
%     incubationTheta(currentTemps>Tm)=0; %reset incubation theta of liquid cells to zero
%     [grid,nuclei,euler,rotMat] = nucleateOctahedrons_v4(grid,nuclei,euler,rotMat,t,currentZHeight,currentTemps,dx,dt,maxTemps,incubationTheta,nuc_params);
%     dTime = toc;
%     timerSaved.findNuclei=timerSaved.findNuclei+dTime;
% 
%     %% Grow the grains in 3D
%     tic
%     [nuclei] = growNuclei_v1(nuclei,currentTemps,Tm,dt,gridSize,realGridSize); %local temperature polynomial growth velocity approach
%     dTime = toc;
%     timerSaved.growNuclei=timerSaved.growNuclei+dTime;
%     %% Save microstructure in vtk format
% %      gridIndicator = grid>0;
% %      binaryGridIndicator = gridIndicator(:) & binaryIndicatory; %find the points that are within the region of interst and also solid
% %      
% %      vtkwrite(sprintf("micro_%i.vtk",iter),'unstructured_grid',xcolumn(binaryGridIndicator),...
% %          ycolumn(binaryGridIndicator),...
% %          zcolumn(binaryGridIndicator),...
% %          'scalars',...
% %          'grain_ID',grid(binaryGridIndicator));

    fprintf("Iteration %i - %.2f%% complete, Simulation %i/%i\n",iter,iter*100/length(beamPath),paramIter,length(params));

    iter=iter+1;
end

%% Print the final structures
iter=iter-1;
xcolumn =xgrid(:);
ycolumn =ygrid(:);
zcolumn =zgrid(:);
tempColumn = currentTemps(:);
binaryIndicatory = zcolumn*realGridSize/gridSize<=beamPath(iter,4); %only print the temperatures at or below the layer height

%  vtkwrite(sprintf("temps_%i.vtk",iter),'unstructured_grid',xcolumn(binaryIndicatory),...
%      ycolumn(binaryIndicatory),...
%      zcolumn(binaryIndicatory),...
%      'scalars',...
%      'temperatures',tempColumn(binaryIndicatory));

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


%% Grain Size Analysis;
%grainDetails = grainSizeAsymmetry(gridSlice,dx,xgrid,ygrid);

%% Print the DAMASK file and save the workspace
writeDAMASK_3D_v1(2,grid,euler,finalGridSize,save_dir); %write the DAMASK files
save(sprintf('workspace_%i.mat',paramIter)); %save the workspace
movefile(sprintf('workspace_%i.mat',paramIter),save_dir+sprintf('/workspace_%i.mat',paramIter)) %mve the file 
movefile(sprintf("micro_%i.vtk",iter),save_dir+sprintf("/micro_%i.vtk",iter)) %mve the file 

end

% %% Plot timers
% figure
% X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Octahedrons','Find New Octahedrons','Remove Octahedrons'});
% Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.findNewOctahedrons,timerSaved.removeNuclei];
% bar(X,Y);
% ylabel("Time [s]");
% 
% %% Plot Pole Figures
% 
cs = crystalSymmetry('cubic');
% rot = rotation.byEuler(euler(:,2),euler(:,3),euler(:,4));
% h = Miller({0,0,1},cs);
% r = rot * h.symmetrise;
% figure
% plot(r)
% title("All grains")
% 
% g = unique(grid(:));
% g=g(g>0);
% 
% rot = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
% r = rot * h.symmetrise;
% figure
% plot(r)
% title("001 - PF, Figure,Original",'contourf')
%mtexColorbar

gridEuler1 = grid(zgrid>80);
g = gridEuler1(:); % = unique(gridEuler1(:));
g=g(g>0);

%rot_mid = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
rot_mid = rotation.byMatrix(rotMat(:,:,g-1000),cs);
h = Miller({0,0,1},cs);
r_mid = rot_mid * h.symmetrise;
figure
plot(r_mid,'antipodal','contourf')
title("001 - PF,Mid Slice",'contourf')
% 
% 
% %% Inverse Pole Figures

g = unique(grid(:));
g=g(g>0);

mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);
figure
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF, Original")

setColorRange([0.6,1.35]) % set equal color range for all subplots
mtexColorbar % add the color bar

figure
plotIPDF(odf,vector3d(0,1,1),'antipodal')
title("011 - IPF, Original")
setColorRange([0.6,1.35]) % set equal color range for all subplots
mtexColorbar % add the color bar

t = norm(odf)^2


% 
% setColorRange([0.6,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% 
% t = norm(odf)^2
% 
% %% Middle slice, should have texture
% 
gridEuler1 = grid(xgrid>50 & xgrid<58);
g = unique(gridEuler1(:));
g=g(g>0);

mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);

figure
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF, Middle")
setColorRange([0.6,2])

t = norm(odf)^2

% setColorRange([0.2,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% figure
% plotIPDF(odf,vector3d(0,1,1),'antipodal')
% title("011 - IPF, Middle")
% 
% setColorRange([0.2,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% t = norm(odf)^2 %texture index - https://mtex-toolbox.github.io/ODFCharacteristics.html
% 
% %% Bottom slice - should be untextured?
% 
% gridEuler1 = grid(zgrid<3);
% g = unique(gridEuler1(:));
% g=g(g>0);
% 
% mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
% odf = unimodalODF(mod1);
% 
% figure
% plotIPDF(odf,vector3d(0,0,1),'antipodal')
% title("001 - IPF, Bottom")
% 
% setColorRange([0.6,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% figure
% plotIPDF(odf,vector3d(0,1,1),'antipodal')
% title("011 - IPF, Bottom")
% 
% setColorRange([0.6,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% t = norm(odf)^2
% 
% 
% %% Random test
% 
randMod1 = rotation.rand(1000);
mod1 = rotation.byEuler(randMod1.phi1,randMod1.Phi,randMod1.phi2,cs);

odf = unimodalODF(mod1);

t = norm(odf)^2
% 
% odf = unimodalODF(mod1);
% figure
% plotIPDF(odf,[zvector],'antipodal')
% title("100 - IPF, Middle")
