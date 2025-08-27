%% Script to model solidification due to a trailing beam path

gridSize = 1000; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

realGridSize=3e-6; %"true" size of the grid. Relevant to the velocity of the boundary motion

poolRadius = 250; %pool radius in grid size units

beamX=500;
beamY=100;

beamdX=0;
beamdY=1;

Temp = 475;
z = 400; % pre-exponential
K1 = 250; % pre-exponential
G_nucleation = 0.25; %Nucleation energy eV/K
Qdiff = 4; %Activation energy for diffusion eV/K
k = 8.617*10^-5; % Boltzmann constant eV/K
Tm=500; %Melting Temp
C = 340; % pre-exponential (temperature independent)
Qgrowth = 0.3; % activation energy for boundary motion/growth eV/K

Gstar = K1*1/(Tm-Temp)^2; %driving force for nucleation due to undercooling
Nhomo=z*exp(-Gstar/(Temp*k))*exp(-G_nucleation/(Temp*k)); %homogeneous nucleation rate

nucleationRate = Nhomo; %per unit time
nucleationRateSTD = 0.1; %nucleation rate STD to use during sampling

growthRate = C*exp(-Qgrowth/(k*Temp));
shapeFactor = 0.005; %heterogeneous nucleation shape factor. If set to anything above 1 then heterogeneous nucleation is suspended

% Timing Constants
totalTime = 1000;
dt = 4;

%Grid to hold the positions
%each x,y position will either be 1 (indicating empty or uncrystallizated space) or a number which indicates the nuclei it is a part of.
grid = ones([gridSize,gridSize]);
grid(1:beamY,gridSize/2-poolRadius:gridSize/2+poolRadius)=-1;

nuclei = []; 
%each row will consist of a unique nuclei
%1: x center position
%2: y center position
%3: time of initilization of nuclei 
%4: color/number assignment of nuclei / grainID
%5: binary indicator of whether the nuclei is still growing or not. 1==growing 0==stopped growing, totally impinged

%Plotting info
cmap = hsv(400); 
cmap = cmap(randperm(400),:);

cmap(1,:)=[1 0 0]; %set the laser color to red
cmap(2,:)=[0 0 0]; %set the liquid color to black
cmap(3,:)=[1 1 1]; %set the 

figure
%Test scanning the beam

[grid,grainCenters] = growEquiaxed(gridSize,200);

%grainCount = length(grainCenters);
%nuclei =[grainCenters(:,1),grainCenters(:,2), zeros([length(grainCenters),1]),(1:length(grainCenters))', zeros([length(grainCenters),1])];
grainCount=1;

for t = 0:dt:totalTime
    %Move the beam
    beamX=beamX+beamdX*dt;
    beamY=beamY+beamdY*dt;
    
    lastGrid=grid;
    if t==0
       %grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius & grid==1)=-1; %move the beam
       grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius)=-1; %move the beam
    else
       %grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius & grid==1)=-1; %move the beam
       grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)<poolRadius)=-1; %move the beam
       
       grid(sqrt((xgrid-beamX).^2+(ygrid-beamY).^2)>poolRadius ...
                & sqrt((xgrid-beamX+beamdX*dt).^2+(ygrid-beamY+beamdY*dt).^2)<poolRadius)=0; %find the newly solidified areas
    end
    
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
    
    %Loop through the growth process for each nuclei
    %Grow all the nuclei
    for i = 1:grainCount-1
        
        %Skip this grain if it is surrounded
        if nuclei(i,5)==0
            continue;
        end
        
        %Find the slice that the grain will grow into
        growthIndices = findNucleiGrowthSlice(xgrid,ygrid,growthRate,nuclei(i,1:2),dt,t,nuclei(i,3));
        
        %Check to see if the grain is surrounded, if so, cease growth
        if ~isempty(growthIndices)
            if isempty(setdiff(growthIndices,find(grid>0)))
                nuclei(i,5)=0;
            end
        end
        
        %If available, assign that region on the grid to the growing grain
        if ~isempty(growthIndices)
            grid(intersect(find(grid==0),growthIndices))=nuclei(i,4);
        end
    end

    %Plot!
    clf
    %imagesc(grid);
    %colormap jet;
    %image(grid,'CDataMapping','direct');
    imshow(grid+2,cmap)
    pause(0.2)
    
end