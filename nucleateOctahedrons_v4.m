function [grid,nuclei,euler,rotMat] = nucleateOctahedrons_v4(grid,nuclei,euler,rotMat,t,topRow,temps,dx,dt,maxTemps,incubationTheta,nuc_params)
%nucleateOctahedrons Function to find the position of new octahedrons
%   ** updated version that newly exposed crystal faces as growing
%   octahedrons
%   ** now including homogeneous nulceation and mushy zone nucleation

Tm = nuc_params.Tm; %melting temperature,K
Nmax = nuc_params.Nmax; %maximum nucleation rate, nuclei / m^3 * s
underCool_N = nuc_params.underCool_N; %mean nucleation undercooing
underCool_sigma = nuc_params.underCool_sigma; %standard deviation nucleation undercooling

Nmax_b = nuc_params.Nmax_b; %maximum nucleation rate for the mushy zone, nuclei / m^3 * s
sigma_b = nuc_params.sigma_b; %standard deviation for mushy zone nucleation, K

%Find all the points that are liquid
liqPoints = find(grid==0);
%currentGridMax = max(max(max(grid))); %highest current index for a grain ID
currentGridMax = max(euler(:,1)); %highest current index for a grain ID

for index = 1:numel(liqPoints)
    
    %Find the location of the index
    [xInd,yInd,zInd] = ind2sub(size(grid),liqPoints(index));

    %Find the neighborhood of points
    [neighborInd] = findNeighbors_3D(xInd,yInd,zInd,grid,topRow);
    neighbors = grid(neighborInd); %find the grain ID of the neighbors

     %If surrounded by liquid, consider the point for homogenous nucleation
    if all(neighbors<1)

         %Find the temperature-dependent nucleation rate
        local_undercooling = Tm-temps(xInd,yInd,zInd);

        %% Bulk Nucleation
        % Temperature-dependent bulk nucleation rate
        P_nucleation = Nmax/(sqrt(2*3.1415)*underCool_sigma)*exp(-((local_undercooling-underCool_N) / (sqrt(2)*underCool_sigma))^2 );
        P_nucleation = P_nucleation*dx^3*dt;

        if rand<P_nucleation %create a new nuclei from bulk homogeneous nucleation
            fprintf("homogeneous zone nucleation\n")
            grainID = currentGridMax+1;
            grid(liqPoints(index))=grainID; %Add the nuclei to the grid
            rot = rotation.rand; %sample a new set of euler angles
            rotMat(:,:,size(rotMat,3)+1) = rot.matrix; %save the rotation matrix
            rotMatrix = rot.matrix;
            rotMatrix=rotMatrix(:);
            rotMatrix=rotMatrix';
            euler = [euler;grainID,rot.phi1,rot.Phi,rot.phi2]; %update the euler index list
            nuclei = [nuclei; %update the nuclei list
            xInd, yInd, zInd,...
            liqPoints(index),...
            rot.phi1,rot.Phi,rot.phi2,...
            t,...
            dx/2,...
            grainID,...
            rotMatrix %rotation matrix
            ];
            currentGridMax = currentGridMax+1;
            continue;
        end
 
        %% Mushy zone nucleation
        %Check if the cell's incubation time is greater than 1
        if incubationTheta(xInd,yInd,zInd)>0

            %Calculate the mushy nucleation probability
            P_nucleation_mush = Nmax_b/(sqrt(2*3.1415)*sigma_b)*exp(-((maxTemps(xInd,yInd,zInd)-Tm) / (sqrt(2)*sigma_b))^2 );
            P_nucleation_mush = P_nucleation_mush*dx^3*dt;

            if rand<P_nucleation_mush %create a new nuclei from mushy zone nucleation
                fprintf("mushy zone nucleation\n")
                grainID = currentGridMax+1;
                grid(liqPoints(index))=grainID; %Add the nuclei to the grid
                rot = rotation.rand; %sample a new set of euler angles
                rotMat(:,:,size(rotMat,3)+1) = rot.matrix; %save the rotation matrix
                rotMatrix = rot.matrix;
                rotMatrix=rotMatrix(:);
                rotMatrix=rotMatrix';                
                euler = [euler;grainID,rot.phi1,rot.Phi,rot.phi2]; %update the euler index list
                nuclei = [nuclei; %update the nuclei list
                xInd, yInd, zInd,...
                liqPoints(index),...
                rot.phi1,rot.Phi,rot.phi2,...
                t,...
                dx/2,...
                grainID,...
                rotMatrix %rotation matrix
                ];
                currentGridMax = currentGridMax+1;
            end     

        end
        continue;

    end

    if sum(neighbors>1)<12 %avoid corner cases of "edge" growth - need to find a better way to manage this?
       continue;
    end

    if size(nuclei,1)>1 %Skip if a nuclei already exists on one of the neighbors  - find a better way to do this?
       if any(neighborInd'==nuclei(:,4))
           continue
       end
    end

    neighbors = neighbors(neighbors>0); %ignore the liquid neighbors

    if isempty(neighbors)
        continue;
    end

    [M,~] = mode(neighbors); %find the mode and frequency
    
    % Epitaxial growth / re-groth of exposed surfaces, 0.00005
%     if rand>0.99 %- setting the minimum to zero activates all of the exposed surfaces... works for IN625
    if rand>0.99 %- setting the minimum to zero activates all of the exposed surfaces...
        %Make sure there aren't multiple modes
        if length(M)==1
                %fprintf("new length(M)==1 nuclei\n")
                eulerIndex = find(euler(:,1)==M);
                grid(liqPoints(index))=M;
                rotationMatrix = rotMat(:,:,eulerIndex);
                rotationMatrix = rotationMatrix(:);
                rotationMatrix = rotationMatrix';
                nuclei = [nuclei;
                    xInd, yInd, zInd,...
                    liqPoints(index),...
                    euler(eulerIndex,2),euler(eulerIndex,3),euler(eulerIndex,4),...
                    t,...
                    0,...
                    M,...
                    rotationMatrix];
        else
                %fprintf("new length(M)>1 nuclei\n")
                eulerIndex = find(euler(:,1)==M(0));
                grid(liqPoints(index))=M(0);
                rotationMatrix = rotMat(:,:,eulerIndex);
                rotationMatrix = rotationMatrix(:);
                rotationMatrix = rotationMatrix';                
                nuclei = [nuclei;
                    XInd, yInd, zInd,...
                    liqPoints(index),...
                    euler(eulerIndex,2),euler(eulerIndex,3),nuclei(euler,4),...
                    t,...
                    0,...
                    M(0),...
                    rotationMatrix];
        end
     end
end


end