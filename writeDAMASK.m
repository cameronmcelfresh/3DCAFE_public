function [] = writeDAMASK(material,grid,xgrid,ygrid,finalGridResolution,dx,cmap,workingDirectory)
%WRITEDAMASK Function to take the results of the AM simulation and 
%write the input files of a DAMASK CP simulation

%% Material Type
%material type 1 == 316L
%material type 2 == Cu
%material type 2> == Aluminum, stand in FCC

if material==1 %for steel
    disDensity = 1e12; 
    G = 79e9; % shear modulus, Pa
    b = 0.254*10^-9; % burgers vector, m
    PN_stress = 200*10^6;
    HP_mult = 20;
elseif material==2 %for Cu
    disDensity = 1e12; 
    G = 88e9; % shear modulus, Pa
    b = 0.26*10^-9; % burgers vector, m
    PN_stress = 15*10^6;
    HP_mult = 3;
else
    disDensity=1e13; % all else
    G = 79e9; % shear modulus, Pa
    b = 0.254*10^9; % burgers vector, m
end


%Make the working directory
mkdir(workingDirectory);

uniqueGrains = unique(grid); %all unique grain IDs

%% Copy the tensionX.yaml file

copyfile('tensionX.yaml',"./"+workingDirectory+"/tensionX.yaml"); %create a fresh file to write over

%% Copy the python geometry prep file

copyfile('geomConvert.py',"./"+workingDirectory+"/geomConvert.py"); %create a fresh file to write over

%% Copy the DAMASK post-processing file for DAMASK

copyfile('DAMASK_postProcess.py',"./"+workingDirectory+"/DAMASK_postProcess.py"); %create a fresh file to write over

%% Construct the material yaml file

%Create a random texture and assign one to each grain 
orients = randrot([length(uniqueGrains),1]);
orients=compact(orients);

%Find the principle axes, ellipse constants, and intersection distance for
%all grains
%figure

pAxes = findprincipleAxes(grid,xgrid,ygrid,[],cmap);

%Calculate the intersection distance for each slip direction on each grain
grid1 = zeros(size(grid));
grid2 = zeros(size(grid));
grid3 = zeros(size(grid));

minDistUniqueGrain = zeros([length(uniqueGrains),1]);

avgDist = 0;
for grainIndex = 1:length(uniqueGrains)
    d_intersect = fit_grainIntersection(quat2rotm(orients(grainIndex,:)),pAxes(grainIndex,4),...
        pAxes(grainIndex,2),pAxes(grainIndex,3),dx);
    
    grid1(grid==uniqueGrains(grainIndex))=d_intersect(1)*dx*10^6; %units in microns!
    grid2(grid==uniqueGrains(grainIndex))=d_intersect(2)*dx*10^6;
    grid3(grid==uniqueGrains(grainIndex))=d_intersect(3)*dx*10^6;

    minDistUniqueGrain(grainIndex) = min(d_intersect)*dx; %save the confinement distance, units in [m]

    %weighted average confinement distance
    avgDist = avgDist+ min(d_intersect)*dx *sum(sum(grid==uniqueGrains(grainIndex)))/(length(grid)*length(grid));
end


% %Plot the results
% figure
% imagesc(grid1)
% title("[1 -1 0] Slip Direction","FontSize",15)
% h = colorbar;
% h.FontSize=15;
% ylabel(h, 'Boundary Confinement Distance [um]')
% caxis([0,60])
% for gI = 1:length(uniqueGrains)
%     hold on
%     quiver([pAxes(gI,9),pAxes(gI,9)],[pAxes(gI,10),pAxes(gI,10)],...
%         [pAxes(gI,5),pAxes(gI,6)],[pAxes(gI,7),pAxes(gI,8)],'k','LineWidth',1);
% end
% axis off
% 
% figure
% imagesc(grid2)
% colorbar
% title("[0 1 1] Slip Direction","FontSize",15)
% h = colorbar;
% h.FontSize=15;
% ylabel(h, 'Boundary Confinement Distance [um]')
% caxis([0,60])
% axis off
% 
% figure
% imagesc(grid3)
% colorbar
% title("[1 0 1] Slip Direction","FontSize",15)
% h = colorbar;
% h.FontSize=15;
% ylabel(h, 'Boundary Confinement Distance [um]')
% caxis([0,60])
% axis off

%%

%Find the grain size of each grain
grainSizes = [];
for grainIndex = 1:length(uniqueGrains)
    %Find the grain size for each grain
    grainArea = dx*dx*sum(sum(grid==uniqueGrains(grainIndex)));
    grainDiameter = sqrt(grainArea/3.1415)*2;
    grainSizes = [grainSizes;grainDiameter];
end

%Output the data necessary to construct the DAMASK inputs 
%writematrix(grainSizes,'grainSizes.txt');
%writematrix(orients,'grainOrientations.txt');

if material==1
    copyfile('material_base_file_316L.txt','material.txt'); %create a fresh file to write over
elseif material==2
     copyfile('material_base_file_Cu.txt','material.txt'); %create a fresh file to write over
else
     copyfile('material_base_file.txt','material.txt'); %create a fresh file to write over
end

fileID = fopen('material.txt','a');

%% Write the yield dependence on grain size

fprintf(fileID,"        xi_0_sl: [%4.3e]\n",G*b*(HP_mult/avgDist+sqrt(disDensity))+PN_stress);
fprintf(fileID,"        output: [xi_sl]\n\n");
fprintf(fileID,"material:\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSERT GRAIN SIZE YIELD DEPENDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(uniqueGrains)
fprintf(fileID,"  - homogenization: SX\n");
fprintf(fileID,"    constituents:\n");

if material==1
    fprintf(fileID,"      - phase: 316L\n");
elseif material==2
     fprintf(fileID,"      - phase: Copper\n");
else
     fprintf(fileID,"      - phase: Aluminum\n");
end
fprintf(fileID,"        v: 1.0\n");
fprintf(fileID,"        xi_0_sl: [%3.5e]\n",G*b*(HP_mult/avgDist+sqrt(disDensity))+PN_stress); % strengthening parameter
%fprintf(fileID,"        gs: %3.4e\n",grainSizes(i));
fprintf(fileID,"        O: [%1.10f, %1.10f,%1.10f, %1.10f]\n",...
    orients(i,1),orients(i,2),orients(i,3),orients(i,4));
end
fclose(fileID);

movefile("material.txt","./"+workingDirectory+"/material.yaml"); %make the film a yaml
%% Write the grainID.geom file

%Downselect the full size of the microstructure
gridDownSelected = imresize(grid,finalGridResolution/length(grid),"nearest");

%relabel the grain ID's to start from 1
gridRelabeled = gridDownSelected;
uniqueRelabeled = unique(gridRelabeled);

for i = 1:length(uniqueRelabeled) %start labeled from 0
    gridRelabeled(gridDownSelected==uniqueRelabeled(i))= i-1;
end

fileID = fopen('grainIDs.txt','w');
fprintf(fileID,'5	header\n');
fprintf(fileID,'grid	a %i	b %i	c 1\n',finalGridResolution, finalGridResolution);
fprintf(fileID,'size	x 1.0	y 1.0	z 0.015625\norigin	x 0.0	y 0.0	z 0.0\nhomogenization	1\n<!skip>\n');
dlmwrite('grainIDs.txt',gridRelabeled,'delimiter',' ','-append');
fclose(fileID);

%Convert grain ID file to the correct postfix (.geom,.yaml)
movefile("grainIDs.txt","./"+workingDirectory+"/grainIDs.geom");

%% Copy the submission scripts
copyfile("submit_simulation.sh","./"+workingDirectory+"/submit_simulation.sh");

end