function tempDist_final = findTempDist3D_v2(rowIter, topRow, sideLength, gridSize, beamPath, tempData)
    %re-written v1 using chatGPT

    % Precompute constants
    heatFluxGridLength = tempData.dx * tempData.gridSize;
    conversionScale = sideLength / heatFluxGridLength;

    beamX = beamPath(rowIter, 2);
    beamY = beamPath(rowIter, 3);

    beamUnitOrientation = beamPath(rowIter, 5:6);
    rotDegree = real(acosd(max(min(dot([1,0], beamUnitOrientation) / (norm([1,0]) * norm(beamUnitOrientation)), 1), -1)));

    tempDist = imrotate_3D(tempData.Temps, rotDegree);

    % Compute indices
    beamX_idx = round(beamX / heatFluxGridLength * tempData.gridSize - 0.5 * conversionScale * tempData.gridSize);
    beamY_idx = round(beamY / heatFluxGridLength * tempData.gridSize - 0.5 * conversionScale * tempData.gridSize);
    
    % Translate beam to true position
    tempDist = circshift(tempDist, [beamY_idx, beamX_idx, 0]);

    % Sample the temp dist matrix
    lengthNeeded = round(tempData.gridSize * conversionScale / 2);
    tempDist_sampled = tempDist(round(tempData.gridSize/2) - lengthNeeded : round(tempData.gridSize/2) + lengthNeeded, ...
                                round(tempData.gridSize/2) - lengthNeeded : round(tempData.gridSize/2) + lengthNeeded, ...
                                end - 2 * lengthNeeded : end);

    % Resize the sampled temperatures to the appropriate length
    tempDist_final = imresizen(tempDist_sampled, gridSize / length(tempDist_sampled));

    % Initialize grid to hold the temperatures
    tempDist_final = zeros(gridSize, gridSize, gridSize);
    % Transfer the melt pool temperatures to the appropriate place in the grid
    tempDist_final(:, :, 1:topRow) = tempDist_sampled(:, :, end - topRow + 1:end);
end
