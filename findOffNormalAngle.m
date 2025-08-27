function rotDegree = findOffNormalAngle(rotMat)
%findOffNormalAngle Function to find the minimum angle of misorientation
%with the central [0,0,1] axis. 

possibleUnitNorms = [0,0,1;
                    0,0,-1;
                    0,1,0;
                    0,-1,0;
                    1,0,0;
                    -1,0,0];

zNormPossible=[0,0,1;
    0,0,-1];

rotDegree=360;
for i = 1:6
    unitNorm = possibleUnitNorms(i,:);
    rotUnit = unitNorm*rotMat; 
    
    for j = 1:2
        zNorm=zNormPossible(j,:);
        CosTheta = max(min(dot(zNorm,rotUnit)/(norm(zNorm)*norm(rotUnit)),1),-1);
        newRotDegree = real(acosd(CosTheta)); 
        if newRotDegree<rotDegree
            rotDegree=newRotDegree;
        end
    end
end

end