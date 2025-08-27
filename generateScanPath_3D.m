function [laserPath] = generateScanPath_3D(velocity, crossHatch, layerThickness, scanRotation, gridSize, timestep, excessLength, startingZ)

%Distance for the laser to step each timestep
deltaX = velocity * timestep;

%Mesh grid of x,y,x points
[Y, X, Z] = meshgrid(-excessLength:crossHatch:(gridSize + excessLength), ...
                     -excessLength:deltaX:(gridSize + excessLength), ...
                     startingZ:layerThickness:(gridSize));        
%Assign the times
time = cumsum(timestep * ones(numel(X), 1),1);

%Placeholder for the unit normals
laserDirectionX = zeros(size(X));
laserDirectionY = zeros(size(X));
laserDirectionX(Z==startingZ)=1;
laserDirectionY(Z==startingZ)=0;

%Unique layer heights
uniqueZ = unique(Z);

%Rotation Matrix
zRotMat = [cosd(scanRotation) -sind(scanRotation) 0;...
            sind(scanRotation) cosd(scanRotation) 0;...
            0   0   1];

%Loop through to apply scan rotations to both the positions and unit normals
for i =2:length(uniqueZ) 
    %Rotate the positions
    xPos = X(Z==uniqueZ(i));
    yPos = Y(Z==uniqueZ(i));

    %Shift center so the positions can be rotated around the central axis
    rotXY = ([xPos(:),yPos(:),zeros(length(xPos),1)] -gridSize/2)*(zRotMat)^(i-1) + gridSize/2;

    if (scanRotation==0 || scanRotation==180) && mod(i,2)==0
        rotXY(:,2)=rotXY(:,2)+crossHatch/2; %add a y-offset shift for every other layer if scan rotation ==0
    end

    X(Z==uniqueZ(i)) = rotXY(:,1);
    Y(Z==uniqueZ(i)) = rotXY(:,2);

    %Rotate the unit vectors for laser direction
    laserDirectionX(Z==uniqueZ(i)) = cosd(-scanRotation*(i-1));
    laserDirectionY(Z==uniqueZ(i)) = sind(-scanRotation*(i-1));       
end

%Vectorize data
laserPath = [time(:), X(:), Y(:), Z(:), laserDirectionX(:), laserDirectionY(:) ];

% % % %Plot, if needed
% figure %laser spots, colored by time
% scatter3(laserPath(:,2),laserPath(:,3),laserPath(:,4),[],laserPath(:,1),'filled');
% xlabel("x")
% ylabel("y")
% 
% figure %quiver plot
% quiver3(laserPath(:,2),laserPath(:,3),laserPath(:,4),laserPath(:,5),laserPath(:,6), zeros(length(laserPath),1));
% xlabel("x")
% ylabel("y")

end

