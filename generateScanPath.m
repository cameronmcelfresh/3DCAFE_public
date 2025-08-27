function [beamPos] = generateScanPath(crossHatch,scanSpeed,realGridSize,dt,beamDiameter,tDelay)
%GENERATESCANPATH Returns the t,x,y position of the beam given a beam path matrix
%with the form of:
%time[s] locationX[% of grid] locationY[% of grid]

% crossHatch == amount of overlap each pass has [%]
% scanSpeed == speed that the beam moves in m/s
% gridSize == size of the grid in pixels
% realGridSize == size of the grid in m
% dt == how often to write out positions for the beam
% delay == amount of time to wait between passes

beamDiameterMultiplier = 1.5;

%beamPos = [0 0 -beamDiameter];
beamPos = [0 0 -beamDiameterMultiplier*beamDiameter];

directVec = 1; %vector to point up and down 
%Calculate the number of ''horizontal'' steps needed and the step distance
stepDist= beamDiameter - beamDiameter*crossHatch;
%stepDist= beamDiameter - beamDiameter*crossHatch*0.5;
stepNum = ceil(realGridSize/stepDist)+1;

for s = 1:stepNum %add a line for each step
    
    %totalDist = realGridSize+beamDiameter*2; %find the necessary vertical distance to travel
    totalDist = realGridSize+beamDiameter*beamDiameterMultiplier*2; %find the necessary vertical distance to travel
    totalTime = totalDist/scanSpeed; %total time to perform the single line scan
    tVec = (max(beamPos(:,1)) :dt : max(beamPos(:,1))+totalTime); %vector of timestamps
        
    initialPos = beamPos(end,2:3); %find the original beam position position
    
    %Check to see if we are above or below the grid
    if initialPos(2)<0
        finalPos = initialPos +[0, totalDist]; %find the final position
    else
        finalPos = initialPos +[0,-totalDist]; %find the final position
    end
    
    inPos = (1:length(tVec))/length(tVec); %fractional position
    inPos_minus = 1-inPos; %fractional position part 2
    
    posVec = initialPos.*inPos_minus' + finalPos.*inPos'; %position vector
    
    %Compile the positions and times
    beamPos = [beamPos;tVec',posVec];
    
    %Add a buffer for the rest of the liquid to solidify
    %Also change the x position of the beam
    beamPos = [beamPos;beamPos(end,1)+tDelay, beamPos(end,2)+stepDist,beamPos(end,3)];
end

end
