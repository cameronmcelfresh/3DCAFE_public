%% Script to model solidification due to a trailing beam path

gridSize = 500; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

realGridSize=3e-6; %"true" size of the grid. Relevant to the velocity of the boundary motion

% Temp = 535;%Temperature in Kelvins
% z = 400; % pre-exponential
% K1 = 250; % pre-exponential
% G_nucleation = 0.15; %Nucleation energy eV/K
% Qdiff = 2.8; %Activation energy for diffusion eV/K
% k = 8.617*10^-5; % Boltzmann constant eV/K
% Tm=500; %Melting Temp
% C = 340; % pre-exponential (temperature independent)
% Qgrowth = 0.28; % activation energy for boundary motion/growth eV/K

Temp = 600;%Temperature in Kelvins
z = 1e14; % pre-exponential
K1 = 3e2; % pre-exponential
G_nucleation = 0.14; %Nucleation energy eV/K
k = 8.617*10^-5; % Boltzmann constant eV/K
C = 7.5e1; % pre-exponential (temperature independent)
Qgrowth = 0.6; % activation energy for boundary motion/growth eV/K

%Tm=1400; %Melting Temp, [K] --> Cu
Tm=1101; %Melting Temp, [K] --> Steel

Gstar = K1*1/(Tm-Temp)^2; %driving force for nucleation due to undercooling
Nhomo=z*exp(-Gstar/(Temp*k))*exp(-G_nucleation/(Temp*k)); %homogeneous nucleation rate

nucleationRate = Nhomo; %per unit time
nucleationRateSTD = 0.1; %nucleation rate STD to use during sampling

growthRate = C*exp(-Qgrowth/(k*Temp));
shapeFactor = 0.8; %heterogeneous nucleation shape factor. If set to anything above 1 then heterogeneous nucleation is suspended

% Timing Constants
timerSaved.initialize=0;
timerSaved.moveBeam=0;
timerSaved.findNuclei=0;
timerSaved.growNuclei=0;
timerSaved.growBoundaries=0;
timerSaved.plot =0;
tic

%Generate the beam path
crossHatch = 0.05; %amount of crossover between scan paths [%]
scanSpeed = 0.3; %speed that the beam rasters [m/s]
beamDiameter = realGridSize/4; % diameter of the beam
poolRadius = round(beamDiameter/(2*realGridSize)*gridSize); %pool radius in grid size units
tDelay = 0.0001; % amount of time to wait between each beam raster
dt = 0.00003; %timestep taken each iteration 

%beam path matrix structure:
%time[s] locationX[% of grid] locationY[% of grid]
beamPath = generateScanPath(crossHatch,scanSpeed,realGridSize,dt,beamDiameter,tDelay);

%Movie Variables
writeMovie=0;
movieTitle="AM_solidification_050_hatch";

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
cmap = hsv(400); 
cmap = cmap(randperm(400),:);

cmap(1,:)=[1 1 1]; %set the laser color to white
cmap(2,:)=[0 0 0]; %set the liquid color to black
cmap(3,:)=[1 1 1]; %set the outside color to white

%Grow the original voronoi tesselation and shift it
[grid,grainCenters] = growEquiaxed(gridSize,200);
grid=grid+100; %shift it so the new grains can start counting at 1

figure

%grainCount = length(grainCenters);
%nuclei =[grainCenters(:,1),grainCenters(:,2), zeros([length(grainCenters),1]),(1:length(grainCenters))', zeros([length(grainCenters),1])];
grainCount=1;
iter=1;

dTime=toc;
timerSaved.initialize = timerSaved.initialize+dTime;

for t = 0:dt:max(beamPath)
    %Move the beam
%     beamX=beamX+beamdX*dt;
%     beamY=beamY+beamdY*dt;
       
      tic
      
      [~,beamInd] = min(abs(beamPath(:,1)-t));
      beamX=beamPath(beamInd,2)/realGridSize*gridSize;
      beamY=beamPath(beamInd,3)/realGridSize*gridSize;
    
    lastGrid=grid;
    if t==0
       %grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius & grid==1)=-1; %move the beam
       grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius)=-1; %move the beam
    else
       grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius)=-1; %move the beam
             
       grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)>poolRadius & ...
           grid<1)=0;
    end
    
    dTime=toc;
    timerSaved.moveBeam = timerSaved.moveBeam+dTime;
    
    tic
    
    %Find all the new nucleation centers
    [centers,centerIndex] = newNucleiPos(grid,xgrid,ygrid,nucleationRate,nucleationRateSTD,dt,shapeFactor);

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
        if nuclei(i,3)==t
            growthIndices = findNucleiGrowthSlice(xgrid,ygrid,growthRate,nuclei(i,1:2),dt,t+dt,nuclei(i,3));
        else
            growthIndices = findNucleiGrowthSlice(xgrid,ygrid,growthRate,nuclei(i,1:2),dt,t,nuclei(i,3));
        end
                        
            
        %Check to see if the grain is surrounded, if so, cease growth
%         if ~isempty(growthIndices)
%             if isempty(setdiff(growthIndices,find(grid>0))) && ~all(grid(growthIndices))==nuclei(i,4)
%                 nuclei(i,5)=0;
%             end
%         end

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
            [growthSlice,growthIndex] = findCAGrowthSlice_timed(xIndex,yIndex,grid,gridSize,filledSpots(i),gridSolidTime, growthRate,realGridSize,t);
            
            %fprintf("%i\n",length(growthSlice));
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
    clf
    imshow(grid+2,cmap)
    pause(0.2)
    if writeMovie 
        %set(gcf,'Position',[15 15 const.movieWidth,const.movieHeight]); %set the position and size of the figure of the figure
        theframe = getframe(gcf);
        videoFrames(iter) = theframe;
    end    
    
    dTime=toc;
    timerSaved.plot=timerSaved.plot+dTime;
    
    iter=iter+1;
end

%% Plot timers

X = categorical({'Initialize','Move Beam','Find Nuclei','Grow Nuclei','Grow Boundaries','Plot'});
Y = [timerSaved.initialize,timerSaved.moveBeam,timerSaved.findNuclei,timerSaved.growNuclei,timerSaved.growBoundaries,timerSaved.plot];

bar(X,Y);
ylabel("Time [s]");


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