%% Script to calculate heat flux of moving source

gridSize = 400;

[xgrid,ygrid] = meshgrid(1:gridSize,1:gridSize);

%Construct beam path
beamPath = (0:1:gridSize)';
beamPath = [beamPath,gridSize/2*ones(size(beamPath))];

dx = 0.005; % step size of the grid [mm]
%heatCap = 377; %heat capacity [J/kg K]
specificHeat = 389; %specific heat 
laserPower = 100; %laser power [W]
v = 200; % scanning velocity [mm/sec]
%interactionTime = 0.05; % interaction time [s]
beamRadius = 0.05; % beam radius [mm]
alpha = 11; % diffusivity [mm^2/s]
density = 8.3600e-06; % density [kg/mm^3]

%The timestep will be set by the velocity and the total distance traveled?
dt = 0.00001; %timestep each iteration [s]

To = 20; %Initial Temperature [C]
Tgrid = To*ones([gridSize,gridSize,length(beamPath)]);
%preFactor = 4*laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

preFactorEager = laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

centerPoint = gridSize/2; %midpoint of the grid
totalSteps = (centerPoint)*dx/(v*dt);

%Find the factor to shift the beam so the center point is the maximum at the end
%of the integration
colShiftFactor = v*dt*totalSteps/dx;

figure
for t=2:totalSteps-2
%     beamX = beamPath(t,1);
%     beamY = beamPath(t,2);

    currentGrid = Tgrid(:,:,t-1);
    for row = 1:gridSize
        for col = 1:gridSize
            
            %Shift the beam so the center point is the maximum at the end
            %of the integration
            col=  col-colShiftFactor/2; %shift the beam for visual purposes
            row = row-gridSize/2;

            %qDot = preFactor*exp( -sqrt((beamX-col).^2+(beamY-row)^2)/(beamRadius^2+8*alpha) );
%             neighbors = findNeighbors(currentGrid,row,col);
%             Tgrid(row,col,t)  = mean(neighbors) + qDot*dx^2/heatCap;

            currentVal = ((totalSteps-t)*dt)^(-0.5)/(beamRadius^2+2*alpha*(totalSteps-t)*dt)*exp( - ( (col*dx-v*t*dt)^2+(row*dx)^2)/(4*alpha *(totalSteps-t)*dt + 2*beamRadius^2) );
            previousVal = ((totalSteps-t-1)*dt)^(-0.5)/(beamRadius^2+2*alpha*(totalSteps-t-1)*dt)*exp( - ( (col*dx-v*(t-1)*dt)^2+(row*dx)^2)/(4*alpha *(totalSteps-t-1)*dt + 2*beamRadius^2));

            if isnan(currentVal)
                fprintf("Found a nan!");
                pause(5);
            end
            
            %qDot = 4*laserPower*exp( -sqrt((beamX-col).^2+(beamY-row)^2)/(beamRadius^2+8*alpha) );
            %preFactorEager = 4*laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

            col=  col+colShiftFactor/2; %shift the beam back to the original location
            row = row+gridSize/2;
            Tgrid(row,col,t) = Tgrid(row,col,t-1) + dt*preFactorEager*(currentVal+previousVal)/2;
            
        end
    end
%     imagesc(Tgrid(:,:,t))
    %contour(Tgrid(:,:,t),5);
%     colorbar
%     caxis([20, 2000]);
%     pause(0.2);
end

%Shift the final matrix
LastGrid = Tgrid(:,:,t);

[maxVals,~]=max(LastGrid);
[~,ind]=max(maxVals); %find maximum indice

LastGrid_shifted= circshift(LastGrid,[0 gridSize/2-ind]);


Tmelt = 1400;
T1 = 1300;
T2 = 1200;
T3 = 1100;
T4 = 1000;
T5 = 900;

figure
imagesc(LastGrid_shifted)
hold on
contour(LastGrid_shifted,[Tmelt,T1,T2,T3,T4,T5],'ShowText','on','Color','w');
colorbar
caxis([20, 2500]);

