function [x,y] = findBeamPos_scanV(t,scanSpeed,beamPath, gridSize)
%findBeamPos Returns the x,y position of the beam given a beam path matrix
%with the form of:
%locationX[% of grid] locationY[% of grid]

%Find the closest time index above the current time

if t>beamPath(end,:) 
    %If the time is beyond the last row of the beam path, set the position
    %to the last listed position
    x = beamPath(end,2)*gridSize;
    y = beamPath(end,3)*gridSize;
else
    minInd = find(beamPath(:,1)-t > 0, 1);

    %Interpolate the nearest x and y position

    %Find the beam positions in fractional coordinates
    xFrac = (beamPath(minInd,1)-t) * beamPath(minInd,2) +  (t-beamPath(minInd-1,1))*beamPath(minInd-1,2);
    yFrac = (beamPath(minInd,1)-t) * beamPath(minInd,3) +  (t-beamPath(minInd-1,1))*beamPath(minInd-1,3);

    %Convert to real coordinates
    x = xFrac*gridSize;
    y = yFrac*gridSize;
end

end

