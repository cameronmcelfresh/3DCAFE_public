function [] = mainSolidification_CA_v1_function(filePath,...
    materialIndicator_array,... %array of indicators for the material type
    scanSpeed_array,... %array of scan speeds to run
    laserPower_array,... %array of laser powers to run
    crossHatch_array,... %array of the crosshatch spacing
    heterogeneousShapeFactor_array, ...%array of shape factors
    duplicates,... %number of duplicates to run
    restartIter) %iteration to restart at
%mainSolidification_CA_v1_function Function to run a solidification trial
%based on input parameters. Output can then be run through DAMASK CPFE for
%mechanical 

%% Material Type
%material type 1 == 316L
%material type 2 == Cu
%material type 2> == Aluminum, stand in FCC

mkdir(filePath)
trialIter=0;

%Loop through all of the combinations
for materialIndicator_Index = 1:length(materialIndicator_array)
for crossHatch_Index = 1:length(crossHatch_array)
for scanSpeed_Index = 1:length(scanSpeed_array)
for laserPower_Index = 1:length(laserPower_array)
for het_Index = 1:length(heterogeneousShapeFactor_array)
for duplicate = 1:duplicates

if trialIter<restartIter
    fprintf("Skipping iteration  " + string(trialIter) + "\n");
    trialIter=trialIter+1;
    continue
end

%Make a folder to hold all of the output data
foldername = filePath+"/trial_"+string(trialIter);
    
%% Assign all of the parametric study variables
materialIndicator = materialIndicator_array(materialIndicator_Index);
crossHatch = crossHatch_array(crossHatch_Index);
scanSpeed = scanSpeed_array(scanSpeed_Index);
laserPower = laserPower_array(laserPower_Index);
het_shapeFactor = heterogeneousShapeFactor_array(het_Index);

%% Run the simulation

gridSize = 150; % smaller grid for test
%gridSize = 500; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid. Relevant to the velocity of the boundary motion

Temp = 600;%Temperature in Kelvins
z = 1e14; % pre-exponential
K1 = 3e2; % pre-exponential
G_nucleation = 0.14; %Nucleation energy eV/K
k = 8.617*10^-5; % Boltzmann constant eV/K
C = 7.5e1; % pre-exponential (temperature independent)
Qgrowth = 0.6; % activation energy for boundary motion/growth eV/K

%Tm=1400; %Melting Temp, [K] --> Cu
Tm=1601; %Melting Temp, [K] --> Steel

Gstar = K1*1/(Tm-Temp)^2; %driving force for nucleation due to undercooling
Nhomo=z*exp(-Gstar/(Temp*k))*exp(-G_nucleation/(Temp*k)); %homogeneous nucleation rate

nucleationRate = Nhomo; %per unit time
nucleationRateSTD = 0.1; %nucleation rate STD to use during sampling

shapeFactor = het_shapeFactor; %heterogeneous nucleation shape factor. If set to anything above 1 then heterogeneous nucleation is suspended

%Save the nucleation info
nucleationData.z = z;
nucleationData.K1 = K1;
nucleationData.G_nucleation = G_nucleation;
nucleationData.k = k;
nucleationData.Tm = Tm;
nucleationData.nucRateSTD = nucleationRateSTD;
nucleationData.shapeFactor=shapeFactor;

growthRate = C*exp(-Qgrowth/(k*Temp));

% Timing Constants
timerSaved.initialize=0;
timerSaved.moveBeam=0;
timerSaved.findNuclei=0;
timerSaved.growNuclei=0;
timerSaved.growBoundaries=0;
timerSaved.plot =0;
tic

%Generate the beam path
%scanSpeed = 0.200; %speed that the beam rasters [m/s]
%laserPower = 100; %laser power in Watts
%crossHatch = 0.1; %amount of crossover between scan paths [%]

%Generate the temperature distribution given the scan speed, laser power,
%and melting temperature
tempData = heat_fluxFunc(scanSpeed*1000,laserPower,Tm);
tempDist = tempData.Temps;

beamDiameter = tempData.poolDiameter; % diameter of the beam
poolRadius = round(beamDiameter/(2*realGridSize)*gridSize); %pool radius in grid size units
tDelay = 0.00001; % amount of time to wait between each beam raster %0.001
%dt = 0.000005; %timestep taken each iteration 
dt = 0.00003; %timestep taken each iteration 


%beam path matrix structure:
%time[s] locationX[% of grid] locationY[% of grid]
beamPath = generateScanPath(crossHatch,scanSpeed,realGridSize,dt,beamDiameter,tDelay);

%Movie Variables
plotData=0; %binary indicator to plot the growing microstructure or not
writeMovie=0;
movieTitle="AM_solidification_test_22_12_27";

if plotData==1
    figure
end

%% Start Computation (all user variables set above)

%Grid to hold the positions
%each x,y position will either be 1 (indicating empty or uncrystallizated space) or a number which indicates the nuclei it is a part of.
grid = ones([gridSize,gridSize]); 
gridSolidTime = zeros([gridSize,gridSize]); %time that the pixel was solidified - for growth of the adjacent grains

nuclei = []; 
%each row will consist of a unique nuclei
%1: x center position
%2: y center position
%3: time of initilization of nuclei 
%4: color/number assignment of nuclei / grainID
%5: binary indicator of whether the nuclei is still growing or not. 1==growing 0==stopped growing, totally impinged

%Pre-compuate all possible boundary points
% indexShaped = reshape(1:numel(grid),[gridSize,gridSize])';
% possibleBoundaries = [indexShaped(beamX-poolRadius,:),indexShaped(beamX+poolRadius,:)];

possibleBoundaries = generateHetBoundaries(grid,beamPath,beamDiameter,realGridSize);

%Plotting info
cmap = hsv(2000); 
cmap = cmap(randperm(2000),:);

cmap(1,:)=[1 1 1]; %set the laser color to white
cmap(2,:)=[0 0 0]; %set the liquid color to black
cmap(3,:)=[1 1 1]; %set the outside color to white

cmapTemp = parula(300);

%Grow the original voronoi tesselation and shift it
[grid,grainCenters] = growEquiaxed(gridSize,200);
grid=grid+1000; %shift it so the new grains can start counting at 1

%grainCount = length(grainCenters);
%nuclei =[grainCenters(:,1),grainCenters(:,2), zeros([length(grainCenters),1]),(1:length(grainCenters))', zeros([length(grainCenters),1])];
grainCount=1;
iter=1;

dTime=toc;
timerSaved.initialize = timerSaved.initialize+dTime;

for t = 0:dt:max(beamPath)
      tic
      
      [~,beamInd] = min(abs(beamPath(:,1)-t));
      currentTemps = findTempDist(beamInd,realGridSize,gridSize,beamPath,tempData);
    
    lastGrid=grid;
    if t==0
       grid(currentTemps>Tm)=-1; %move the beam
    else
       grid(currentTemps>Tm)=-1; %move the beam
             
       grid(currentTemps<Tm & ...
           grid<1)=0;
    end
    
    dTime=toc;
    timerSaved.moveBeam = timerSaved.moveBeam+dTime;
    
    tic
    
    %Find all the new nucleation centers
    %[centers,centerIndex] = newNucleiPos(grid,xgrid,ygrid,nucleationRate,nucleationRateSTD,dt,shapeFactor);
    [centers,centerIndex] = newNucleiPos_CA(grid,xgrid,ygrid,currentTemps,nucleationData,dt,realGridSize/gridSize);

    grainCountOrig=grainCount; %save the grain count variable to be reused in the next loop
    
    %Append the new nuclei to the nuclei matrix
    for i = 1:length(centerIndex)
        grainCount=grainCount+1;
        nuclei=[nuclei;...
            centers(i,1),centers(i,2),t,grainCount,1];
    end
    
    %Add the nuclei to the grid
    for i=1:length(centerIndex)
        grainCountOrig=grainCountOrig+1;
        grid(centerIndex(i))=grainCountOrig;
    end    
    
    dTime=toc;
    timerSaved.findNuclei=timerSaved.findNuclei+dTime;
    
    tic
    
    %Loop through the growth process for each nuclei
    %Grow all the nuclei
    for i = 1:grainCount-1
        
        %Skip this grain if it is surrounded
        if nuclei(i,5)==0
            continue;
        end
        
        %Find the slice that the grain will grow into          
        localTempData = currentTemps(grid==nuclei(i,4)); %Find the local temperatures
        localTemp = mean(localTempData(:)); %find the average temperature
        growthRateLocal = C*exp(-Qgrowth/(k*localTemp));
        
        %Ensure the local growth rate is below the scan speed
        if growthRateLocal>scanSpeed*0.6
            growthRateLocal = scanSpeed*0.6;
        elseif growthRateLocal<scanSpeed*0.1
            growthRateLocal = scanSpeed*0.1;
        end

        %If the nuclei doesn't have any points on the grid, remove it
        if sum(sum(grid==nuclei(i,4)))==0
            nuclei(i,5)=0;
            continue
        end

        growthIndices = findNucleiGrowthSlice_1(xgrid,ygrid,growthRateLocal,nuclei(i,1:2),dt,realGridSize/gridSize,grid,nuclei(i,4));            

        if isSurrounded(i+1,grid) %the plus 1 is needed because the grainID's start counting from 2
            nuclei(i,5)=0;
        end
        
        %If available, assign that region on the grid to the growing grain
        if ~isempty(growthIndices)
            grid(intersect(find(grid==0),growthIndices))=nuclei(i,4);
        end
    end
    
    dTime=toc;
    timerSaved.growNuclei=timerSaved.growNuclei+dTime;
    
    tic
    
    %Find all the boundaries
    filledSpots = intersect(find(grid>2),possibleBoundaries);
    
    %Find possible boundaries - use cellular automata growth
    for i=1:length(filledSpots)
        
        xIndex = xgrid(filledSpots(i));
        yIndex = ygrid(filledSpots(i));
        
        %Update the surrounding spots if it is boundary pixel
        if isBoundary(xIndex,yIndex,grid) 
            
            %Initiate the boundary wtih a solidification time
            if gridSolidTime(filledSpots(i))==0
                gridSolidTime(filledSpots(i))=t;
            end
            
            %Find the slice that the area would grow into
            %[growthSlice,growthIndex] = findCAGrowthSlice(xIndex,yIndex,grid,gridSize,filledSpots(i)); %non-physical way
            %[growthSlice,growthIndex] =
            %findCAGrowthSlice_timed(xIndex,yIndex,grid,gridSize,filledSpots(i),gridSolidTime,
            %growthRate,realGridSize,t); %physical method with constant
            %background temperature
            growthRateLocal = C*exp(-Qgrowth/(k*currentTemps(yIndex,xIndex)))*gridSize/realGridSize; %local growth rate in units of cell step /s
            
            %Ensure the local growth rate is below the scan speed
            if growthRateLocal>scanSpeed*gridSize/realGridSize*0.6
                growthRateLocal = scanSpeed*0.6*gridSize/realGridSize;
            elseif growthRateLocal<scanSpeed*gridSize/realGridSize*0.1
                growthRateLocal = scanSpeed*0.1*gridSize/realGridSize;
            end
            [growthSlice,growthIndex] = findCAGrowthSlice_timed(xIndex,yIndex,grid,gridSize,filledSpots(i),gridSolidTime, growthRateLocal,realGridSize,t);
            
            %Assign the points on the grid
            grid(growthIndex)=grid(filledSpots(i));
            
            %Save the solidifying time
            gridSolidTime(growthIndex)=t;
            
            %Increment the voxels of boundaries
            possibleBoundaries=[possibleBoundaries,growthIndex'];    
        end
    end
    
    dTime=toc;
    timerSaved.growBoundaries=timerSaved.growBoundaries+dTime;

    tic 
    
    %Plot!
    if plotData==1
        subplot(1,2,1); %plot the grain structure
        imshow(grid+2,cmap)
        axis off
        subplot(1,2,2); %plot the laser temperature profile
        imagesc(currentTemps)
        colorbar
        caxis([20, 2000]);
        axis off
        set(gcf, 'Position',  [100, 100, 1300, 550])
        pause(0.0001)
        if writeMovie 
            %set(gcf,'Position',[15 15 const.movieWidth,const.movieHeight]); %set the position and size of the figure of the figure
            theframe = getframe(gcf);
            videoFrames(iter) = theframe;
        end    
    end
    
    dTime=toc;
    timerSaved.plot=timerSaved.plot+dTime;
    
    iter=iter+1;
end

%% Plot timers

if plotData==1
    X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Nuclei','Grow Boundaries','Plot'});
    Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.growBoundaries,timerSaved.plot];
    bar(X,Y);
    ylabel("Time [s]");
end

%% Save the video
if writeMovie==1 %Record a movie if specified to do so   
    moviename=sprintf(movieTitle+'.avi'); %update the name of the file
    aviobj=VideoWriter(moviename);
    aviobj.Quality=100;

    open(aviobj);    
    for frameIter = 1:length(videoFrames)
        writeVideo(aviobj,videoFrames(frameIter));
    end
    close(aviobj);
end

%% Smooth
grid = modeFilter(grid,4);

%% Write the DAMASk Files
writeDAMASK(materialIndicator,grid,xgrid,ygrid,200,realGridSize/gridSize,cmap,foldername)

%% Save the workspace
mkdir(foldername);
filename = 'mat_dat.mat';
save(foldername + "/" + filename);

%% Save relevant parametric study data
fid = fopen( foldername + '/parametricDetails.txt', 'wt' );

materialIndicator = materialIndicator_array(materialIndicator_Index);
crossHatch = crossHatch_array(crossHatch_Index);
scanSpeed = scanSpeed_array(scanSpeed_Index);
laserPower = laserPower_array(laserPower_Index);
het_shapeFactor = heterogeneousShapeFactor_array(het_Index);

fprintf( fid, 'materialIndicator, %f\n', materialIndicator);
fprintf( fid, 'crossHatch, %f\n', crossHatch);
fprintf( fid, 'scanSpeed, %f\n', scanSpeed);
fprintf( fid, 'laserPower, %f\n', laserPower);
fprintf( fid, 'het_shapeFactor, %f\n', het_shapeFactor);

fclose(fid);

trialIter=trialIter+1; %add to the iterator and move to the next trial

end %end duplicates
end %end hetergeneous index
end %end laser power index
end %end scan speed index
end %end cross hatch index
end %end material index

%% Write the group submission file
fid = fopen( filePath + '/group_submit.sh', 'wt' );

for iter_val = 1:trialIter-1
    fprintf( fid, 'cd trial_%i\n', iter_val);
    fprintf( fid, './submit_simulation.sh\ncd ..\n\n');
end

end