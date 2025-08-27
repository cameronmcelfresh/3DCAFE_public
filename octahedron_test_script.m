%Script to test the decentered octahedron approach

gridSize = 100; % side length for square grid of size (gridSize,gridSize) when first constructing the vertices
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));

realGridSize=300e-6; %"true" size of the grid [m]. Relevant to the velocity of the boundary motion

aOctahedron = zeros([gridSize gridSize gridSize]);

euler_angles=[0,0,0];

grid(50,50,50)=1;
center = [50 50 50];
solidTime = zeros([grid, grid,grid]);

growthRate = 0.001;
dt = 0.0001;
totalTime = dt*100;

octaHedrons = [50,50,50,0,4950];
%x,y,z,t,i
%x/y/z position
%t=time solidified
%i=index of center point for mushy cell

for t=0:dt:totalTime

    for x = 1:gridSize
        for y = 1:gridSize
            for z = 1:gridSize
    
            if is_inside_rotated_octahedron(x, y, z, center,side_length, euler_angles)
                aOctahedron(x,y,z)=1;
            end
    
            end
        end
    end
end
binaryIndicatory = aOctahedron==1;

binaryIndicatory=binaryIndicatory(:);
aOctahedron=aOctahedron(:);

xgrid=xgrid(:);
ygrid=ygrid(:);
zgrid=zgrid(:);

vtkwrite("octahedron.vtk",'structured_grid',xgrid(binaryIndicatory),...
    ygrid(binaryIndicatory),...
    zgrid(binaryIndicatory),...
    'scalars',...
    'ID',aOctahedron(binaryIndicatory));