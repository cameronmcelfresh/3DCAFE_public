%% Scan path test

%% Script to model solidification due to a trailing beam path

gridSize = 1000; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

realGridSize=800e-6; %"true" size of the grid. Relevant to the velocity of the boundary motion

Temp = 535;%Temperature in Kelvins
Tm=1600; %Melting Temperature in Kelvin
scanSpeed = 0.100; %speed that the beam rasters [m/s]
laserPower = 200; %laser power in Watts

%Generate the temperature distribution given the scan speed, laser power,
%and melting temperature
tempData = heat_fluxFunc(scanSpeed*1000,laserPower,Tm);
tempDist = tempData.Temps;


%Generate the beam path
crossHatch = 0; %amount of crossover between scan paths [%]
%beamDiameter = realGridSize/4; % diameter of the beam
beamDiameter = tempData.poolDiameter; % diameter of the beam
poolRadius = round(beamDiameter/(2*realGridSize)*gridSize); %pool radius in grid size units
tDelay = 150; % amount of time to wait between each beam raster
dt = 0.0001; %timestep taken each iteration


%beam path matrix structure:
%time[s] locationX[% of grid] locationY[% of grid]
beamPath = generateScanPath(crossHatch,scanSpeed,realGridSize,dt,beamDiameter,tDelay);

%Movie Variables
writeMovie=0;
movieTitle="ET_p200_v50";

%% Start Computation (all user variables set above)

%Grid to hold the positions
%each x,y position will either be 1 (indicating empty or uncrystallizated space) or a number which indicates the nuclei it is a part of.
grid = ones([gridSize,gridSize]); 
gridSolidTime = zeros([gridSize,gridSize]); %time that the pixel was solidified - for growth of the adjacent grains

figure

iter=1;

%multFactor=200;

%for t = 0:dt:max(beamPath)
for t = 0:dt:max(beamPath(:,1))
    [~,beamInd] = min(abs(beamPath(:,1)-t));
    %beamX=round(beamPath(beamInd,2)/realGridSize*gridSize);
    %beamY=round(beamPath(beamInd,3)/realGridSize*gridSize);
    %beamX=beamPath(beamInd,2);
    %beamY=beamPath(beamInd,3);
        
    currentTemps = findTempDist(iter,realGridSize,gridSize,beamPath,tempData);
    
    %Plot!
    clf
    imagesc(currentTemps)
    colorbar
    caxis([20, 2000]);
    axis off
    pause(0.01);
    if writeMovie 
        %set(gcf,'Position',[15 15 const.movieWidth,const.movieHeight]); %set the position and size of the figure of the figure
        theframe = getframe(gcf);
        videoFrames(iter) = theframe;
    end    

    iter=iter+1;
end


%% Save the video
if writeMovie==1 %Record a movie if specified to do so   
    moviename=sprintf(movieTitle+'.avi'); %update the name of the file
    aviobj=VideoWriter(moviename);
    aviobj.Quality=100;

    open(aviobj);    
    for frameIter = 1:length(videoFrames)
        try
            writeVideo(aviobj,videoFrames(frameIter));
        end
    end
    close(aviobj);
end