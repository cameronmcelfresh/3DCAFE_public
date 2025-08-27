%Script to test the decentered octahedron approach

gridSize = 100; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion

grid = zeros([gridSize gridSize gridSize]);

grid(40,50,50)=1;
grid(60,50,50)=2;

growthRate = 5000;
dt = 0.0001;
totalTime = dt*60;

rotMat = eye(3);
rot = rotation.rand(1);
rotMat1 = rot.matrix; %rotation matrices to handle the octahedron rotations. Due to MTEX Bunge rotations, the transpose of the rotation matrix is needed for proper euler configuration.

%nuclei = [50,50,50,494950,0,0,0,0,0,1];
nuclei = [40,50,50,494940,0,0,0,0,0,1,rotMat(:)';...
            60,50,50,494960,rot.phi1,rot.Phi,rot.phi2,0,0,2,rotMat1(:)'];

% column 1-3: x,y,z position of the nuclei
% column 4: index of the cell that the nuclei belongs to
% column 5-7: euler angles of the nuclei
% column 8: initial nucleation time of the nuclei
% column 9: size of the nuclei
% column 10: grainID of nuclei

iter=1;

for t=0:dt:totalTime

    %Remove any surrounded nuclei
    nuclei = removeSurroundedNuclei(nuclei,grid);

    %Find new octohedrons
    %[grid,nuclei] = findNewOctohedrons(grid,nuclei,t);
    [grid,nuclei] = findNewOctohedrons_v2(grid,nuclei,t,10000);

    %Print proress
    binaryIndicatory = grid>0;

    binaryIndicatory=binaryIndicatory(:);
    gridColumn=grid(:);
    
    xgridColumn=xgrid(:);
    ygridColumn=ygrid(:);
    zgridColumn=zgrid(:);
    
    vtkwrite(sprintf("octahedron_%i.vtk",iter),'unstructured_grid',xgridColumn(binaryIndicatory),...
        ygridColumn(binaryIndicatory),...
        zgridColumn(binaryIndicatory),...
        'scalars',...
        'grainID',gridColumn(binaryIndicatory));

    %Grow nuclei
    [nuclei] = growNuclei(nuclei,growthRate,dt);

    iter=iter+1;
end

