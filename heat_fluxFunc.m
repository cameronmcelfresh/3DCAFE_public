function [tempData] = heat_fluxFunc(v,laserPower,Tm)
%Function to calculate the temperature profile and melt pool dimensions
%Tm == melting temperature in Kelvin

%Structure to hold the data
tempData = {};

gridSize = 1000;

%Construct beam path
beamPath = (0:1:gridSize)';
beamPath = [beamPath,gridSize/2*ones(size(beamPath))];

dx = 0.001; % step size of the grid [mm]
%beamRadius = 0.025; % beam radius [mm]

% Constants for Cu
% specificHeat = 389; %specific heat --> Cu
% alpha = 50; % diffusivity [mm^2/s] --> fit from Cu laser paper
% density = 8.3600e-06; % density [kg/mm^3] --> Cu
% beamRadius = 0.025; % beam radius [mm] --> Cu

% Constants for 316L
specificHeat = 770; %specific heat --> 316L
alpha = 20; % diffusivity [mm^2/s] --> 316L
density = 8.0000e-06; % density [kg/mm^3] --> 316L
%beamRadius = 0.005; % beam radius [mm] --> cont. used for 2023_04_07_parametric study steel
beamRadius = 0.025; % beam radius [mm] --> 

%The timestep will be set by the velocity and the total distance traveled?
dt = 0.00001; %timestep each iteration [s]

To = 20; %Initial Temperature [C]
Tgrid = To*ones([gridSize,gridSize]);
%Tgrid = To*ones([gridSize,gridSize,length(beamPath)]);
%preFactor = 4*laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

preFactorEager = laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

centerPoint = gridSize/2; %midpoint of the grid
totalSteps = (centerPoint)*dx/(v*dt);

%Find the factor to shift the beam so the center point is the maximum at the end
%of the integration
colShiftFactor = v*dt*totalSteps/dx;

for t=2:totalSteps-2
%     beamX = beamPath(t,1);
%     beamY = beamPath(t,2);

    %currentGrid = Tgrid(:,:,t-1);
    for row = 1:gridSize
        for col = 1:gridSize
            
            %Shift the beam so the center point is the maximum at the end
            %of the integration
            col=  col-colShiftFactor/2; %shift the beam for visual purposes
            row = row-gridSize/2;

            currentVal = ((totalSteps-t)*dt)^(-0.5)/(beamRadius^2+2*alpha*(totalSteps-t)*dt)*exp( - ( (col*dx-v*t*dt)^2+(row*dx)^2)/(4*alpha *(totalSteps-t)*dt + 2*beamRadius^2) );
            previousVal = ((totalSteps-t-1)*dt)^(-0.5)/(beamRadius^2+2*alpha*(totalSteps-t-1)*dt)*exp( - ( (col*dx-v*(t-1)*dt)^2+(row*dx)^2)/(4*alpha *(totalSteps-t-1)*dt + 2*beamRadius^2));

            if isnan(currentVal)
                fprintf("Found a nan!");
                pause(5);
            end
            
            col=  col+colShiftFactor/2; %shift the beam back to the original location
            row = row+gridSize/2;
            Tgrid(row,col) = Tgrid(row,col) + dt*preFactorEager*(currentVal+previousVal)/2;
            %Tgrid(row,col,t) = Tgrid(row,col,t-1) + dt*preFactorEager*(currentVal+previousVal)/2;
            
        end
    end
%     imagesc(Tgrid(:,:,t))
    %contour(Tgrid(:,:,t),5);
%     colorbar
%     caxis([20, 2000]);
%     pause(0.2);
end

%Shift the final matrix
%LastGrid = Tgrid(:,:,t);
LastGrid = Tgrid;

[maxVals,~]=max(LastGrid);
[~,ind]=max(maxVals); %find maximum indice

%shift the temperature grid to the center point
LastGrid_shifted= circshift(LastGrid,[0 gridSize/2-ind]); 

%Calculate the pool width and diameter
%Length
poolLength = dx/1000*sum(max(LastGrid_shifted,[],1)>Tm);
poolDiameter = dx/1000*sum(max(LastGrid_shifted,[],2)>Tm);
 
tempData.Temps = LastGrid_shifted; %temperatures
tempData.dx = dx/1000; %mm -> m conversion
tempData.gridSize =gridSize; %grid size
tempData.To = To; %base temperature
tempData.poolLength = poolLength; %melt pool length in m
tempData.poolDiameter = poolDiameter; %melt pool diameter in m

end

