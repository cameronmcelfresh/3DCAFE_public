%% Script to model solidification due to a trailing beam path
startup_mtex; %start mtex

save_dir = "test";
finalGridSize = [50,50]; %Grid size to simulate with CPFE
%%
gridSize = 125; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion
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
crossHatch = 85e-6; %amount to step laterally each scan path [m]
layerThickness = 40e-6; %layer thickness during AM [m]
excessLength = 110e-6; %distance to start outside of the grid [m]
scanRotation=0; % degrees to rotate after each rotation [degrees]
startingZ = 80e-6; %height to start the first layer [m]
numGrains = 40; %number of grains to initialize in the voronoi tesselation

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
tDelay = 0.000005; % amount of time to wait between each beam raster %0.000025
dt = 0.0000015; %timestep taken each iteration , dt = 0.00005;

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
nuc_params.Tm = 1400; %melting temperature,K
nuc_params.Nmax = 5e20; %maximum nucleation rate, nuclei / m^3 * s
nuc_params.underCool_N = 5; %mean nucleation undercooing
nuc_params.underCool_sigma = 1; %standard deviation nucleation undercooling

nuc_params.Nmax_b = 5e22; %maximum nucleation rate for the mushy zone, nuclei / m^3 * s
nuc_params.sigma_b = 50; %standard deviation for mushy zone nucleation, K

Cinc = 2e-3; %incubation time constant, K*s
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
    [grid,nuclei] = findNewOctohedrons_v2(grid,nuclei,t,currentZHeight);
    dTime = toc;
    timerSaved.findNewOctahedrons=timerSaved.findNewOctahedrons+dTime;

    %% Nucleate new octahedrons
    tic
    incubationTheta = incubationTheta + dt./(Cinc./(Tm-currentTemps)); %update the incubation time matrix
    incubationTheta(grid>0)=0; %reset incubation theta of solid cells to zero
    incubationTheta(currentTemps>Tm)=0; %reset incubation theta of liquid cells to zero
    [grid,nuclei,euler,rotMat] = nucleateOctahedrons_v4(grid,nuclei,euler,rotMat,t,currentZHeight,currentTemps,dx,dt,maxTemps,incubationTheta,nuc_params);
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

    fprintf("Iteration %i - %.2f%% complete\n",iter,iter*100/length(beamPath));

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


%% Grain Size Analysis;
%grainDetails = grainSizeAsymmetry(gridSlice,dx,xgrid,ygrid);

%% Print the DAMASK file and save the workspace
writeDAMASK_3D_v1(1,grid,euler,finalGridSize,save_dir); %write the DAMASK files
save(save_dir+'\workspace.mat'); %save the workspace
movefile(sprintf("micro_%i.vtk",iter),save_dir+sprintf("/micro_%i.vtk",iter)) %mve the file 

%%

% %% Plot timers
% figure
% X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Octahedrons','Find New Octahedrons','Remove Octahedrons'});
% Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.findNewOctahedrons,timerSaved.removeNuclei];
% bar(X,Y);
% ylabel("Time [s]");
% 
% %% Plot Pole Figures
% 
% cs = crystalSymmetry('cubic');
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
% %mtexColorbar
% 
% gridEuler1 = grid(xgrid>40 & xgrid<75);
% g = unique(gridEuler1(:));
% g=g(g>0);
% 
% rot_mid = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
% 
% r_mid = rot_mid * h.symmetrise;
% figure
% plot(r_mid)
% title("001 - PF,Mid Slice",'contourf')
% 
% 
% %% Inverse Pole Figures
% 
% g = unique(grid(:));
% g=g(g>0);
% 
% mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
% odf = unimodalODF(mod1);
% figure
% plotIPDF(odf,vector3d(0,0,1),'antipodal')
% title("001 - IPF, Original")
% 
% setColorRange([0.6,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% figure
% plotIPDF(odf,vector3d(0,1,1),'antipodal')
% title("011 - IPF, Original")
% 
% setColorRange([0.6,2]) % set equal color range for all subplots
% mtexColorbar % add the color bar
% 
% 
% t = norm(odf)^2
% 
% %% Middle slice, should have texture
% 
% gridEuler1 = grid(zgrid>80 & zgrid<82);
% g = unique(gridEuler1(:));
% g=g(g>0);
% 
% mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
% odf = unimodalODF(mod1);
% 
% figure
% plotIPDF(odf,vector3d(0,0,1),'antipodal')
% title("001 - IPF, Middle")
% 
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
% randMod1 = rotation.rand(1000);
% mod1 = rotation.byEuler(randMod1.phi1,randMod1.Phi,randMod1.phi2,cs);
% 
% odf = unimodalODF(mod1);
% 
% t = norm(odf)^2
% 
% odf = unimodalODF(mod1);
% figure
% plotIPDF(odf,[zvector],'antipodal')
% title("100 - IPF, Middle")
