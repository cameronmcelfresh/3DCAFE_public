function [nuclei] = growNuclei_v1(nuclei,temps,Tm,dt,gridSize,realGridSize)
%growNuclei Function to grow all of the nuclei
%  updated function that includes local temperature effects
%  growth rates are calculated from the local undercooling conditions

% column 1-3: x,y,z position of the nuclei
% column 4: index of the cell that the nuclei belongs to
% column 5-7: euler angles of the nuclei
% column 8: initial nucleation time of the nuclei
% column 9: size of the nuclei


%% Growth constants for IN625

A = -1.0302*10^-7; 
B = 1.0533*10^-4;
C = 2.2196*10^-3; 
D = 0;

%% Growth constants for 316L

% A = 1.091*10^-5;
% B = -2.034*10^-4;
% C = 2.74*10^-3; 
% D = 1.151*10^-4;

A = A/2;
B = B/2;
C = C/2; 
D = D/2;

%% Perform growth for each 

if size(nuclei,1)>0 %only perform if there are any nuclei to grow
    
    %Find the cells that each nuclei belongs to. Enforce bounds if needed
    nucleiPos = [round(nuclei(:,1)),round(nuclei(:,2)),round(nuclei(:,3))];
    nucleiPos(nucleiPos<1) = 1;
    nucleiPos(nucleiPos>length(temps)) = length(temps);

    %Find the local temperature in K
    localInd = sub2ind(size(temps), nucleiPos(:,1), nucleiPos(:,2), nucleiPos(:,3)); %rough index of octahedrons
    localT = temps(localInd); %temperature of the cell
    underCooling = Tm-localT; %undercooling of the cell

    %Calculate the local growth rate
    growthV = A*underCooling.^3 + B.*underCooling.^2 + C*underCooling + D;
    growthV = growthV * gridSize/realGridSize; % Convert from "true" growth rate to "pixel" growth rate

    growthV(growthV<0) = 0; % Only grow the nuclei if the growth rate is positive
    growthV(underCooling<0) = 0; % Avoid growth at high temperatures

    %nuclei(:,9)= nuclei(:,9) + growthV*dt;
    nuclei(:,9)= nuclei(:,9) + growthV*dt/10;
end

end