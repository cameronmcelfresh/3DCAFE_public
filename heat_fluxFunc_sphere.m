function [tempData] = heat_fluxFunc_sphere(v,laserPower,Tm,dx)
%Function to calculate the temperature profile and melt pool dimensions
%Tm == melting temperature in Kelvin

%Structure to hold the data
tempData = {};

gridSize = 400;

[xgrid,ygrid,zgrid]=meshgrid(1:gridSize,1:gridSize,1:gridSize);

xgrid=xgrid(:);
ygrid=ygrid(:);
zgrid=zgrid(:);

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
Tgrid = To*ones([gridSize,gridSize,gridSize]);
%Tgrid = To*ones([gridSize,gridSize,length(beamPath)]);
%preFactor = 4*laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

preFactorEager = laserPower/(density*specificHeat*3.14*sqrt(4*alpha*3.14));

centerPoint = gridSize/2; %midpoint of the grid
totalSteps = 2;

%Find the factor to shift the beam so the center point is the maximum at the end
%of the integration

for t=1:totalSteps

    %currentGrid = Tgrid(:,:,t-1);
    for row = 1:gridSize
        for col = 1:gridSize
            for depth =1: gridSize

                dist = sqrt( (row-gridSize/2)^2 + (col-gridSize/2)^2 + (depth-gridSize)^2)*1e-6;
                T=300;
                if dist<50e-6
                    T=1500;
                end

                Tgrid(row,col,depth) = T;
            end
        end
    end

end

%Tunravel = Tgrid(:);
%scatter3(xgrid(1:13:end),ygrid(1:13:end),zgrid(1:13:end),[],Tunravel(1:13:end))
%colorbar
%clim([300, 1500]);

%Shift the final matrix
LastGrid = Tgrid;

%[maxVals,~]=max(LastGrid);
%[~,ind]=max(maxVals); %find maximum indice

%shift the temperature grid to the center point
%LastGrid_shifted= circshift(LastGrid,[0 gridSize/2-ind gridSize/2-ind]); 
%LastGrid_shifted= circshift(LastGrid,[0 gridSize/2-ind gridSize/2-ind]); 
LastGrid_shifted= LastGrid;


%Calculate the pool width and diameter
%Length
poolLength = 50e-6 ; %dx/1000*sum(max(LastGrid_shifted,[],1)>Tm);
poolDiameter = 50e-6; %dx/1000*sum(max(LastGrid_shifted,[],2)>Tm);
 
tempData.Temps = LastGrid_shifted; %temperatures
tempData.dx = dx/1000; %mm -> m conversion
tempData.gridSize =gridSize; %grid size
tempData.To = To; %base temperature
tempData.poolLength = poolLength; %melt pool length in m
tempData.poolDiameter = poolDiameter; %melt pool diameter in m

end

