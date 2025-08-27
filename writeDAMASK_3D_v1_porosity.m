function [gridPorous,gridDownSelected1] = writeDAMASK_3D_v1_porosity(material,grid,xgrid,ygrid,zgrid,Tm,euler,finalGridResolution,maxTempMat,workingDirectory)
%WRITEDAMASK Function to take the results of the AM simulation and 
%write the input files of a DAMASK CP simulation - using porosity

%% Material Type
%material type 1 == 316L
%material type 2 == IN625
%material type 2> == Aluminum, stand in FCC

if material==1 %for steel
    disDensity = 1e12; 
    G = 79e9; % shear modulus, Pa
    b = 0.254*10^-9; % burgers vector, m
    PN_stress = 200*10^6;
    HP_mult = 20;
elseif material==2 %for IN625
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

%%

gridPorous = createPorosity(grid,xgrid,ygrid,zgrid,Tm*1.1,maxTempMat,30);

%%

%Make the working directory
mkdir(workingDirectory);

uniqueGrains = unique(gridPorous); %all unique grain IDs

%% Copy the tensionX.yaml file

copyfile('tensionZ.yaml',"./"+workingDirectory+"/tensionZ.yaml"); %create a fresh file to write over
copyfile('tensionX.yaml',"./"+workingDirectory+"/tensionX.yaml"); %create a fresh file to write over
copyfile('tensionY.yaml',"./"+workingDirectory+"/tensionY.yaml"); %create a fresh file to write over

%% Copy the python geometry prep file

copyfile('geomConvert.py',"./"+workingDirectory+"/geomConvert.py"); %create a fresh file to write over

%% Copy the DAMASK post-processing file for DAMASK

%copyfile('DAMASK_postProcess.py',"./"+workingDirectory+"/DAMASK_postProcess.py"); %create a fresh file to write over
copyfile('postProcessXX.py',"./"+workingDirectory+"/postProcessXX.py"); %create a fresh file to write over
copyfile('postProcessYY.py',"./"+workingDirectory+"/postProcessYY.py"); %create a fresh file to write over
copyfile('postProcessZZ.py',"./"+workingDirectory+"/postProcessZZ.py"); %create a fresh file to write over

%% Construct the material yaml file

%%
%Output the data necessary to construct the DAMASK inputs 

if material==1
    copyfile('material_base_file_316L_voids.txt','material.txt'); %create a fresh file to write over
elseif material==2
     copyfile('material_base_file_IN625_voids.txt','material.txt'); %create a fresh file to write over
else
     copyfile('material_base_file.txt','material.txt'); %create a fresh file to write over
end

fileID = fopen('material.txt','a');

%% Write the yield dependence on grain size
%avgDist = 10e-6;
%
%fprintf(fileID,"        xi_0_sl: [%4.3e]\n",G*b*(HP_mult/avgDist+sqrt(disDensity))+PN_stress);
%fprintf(fileID,"        output: [xi_sl]\n\n");
%fprintf(fileID,"material:\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSERT GRAIN SIZE YIELD DEPENDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Print the void properties
fprintf(fileID,"  - homogenization: SX\n");
fprintf(fileID,"    constituents:\n");
fprintf(fileID,"      - phase: void\n");
fprintf(fileID,"        v: 1.0\n");
quat = eul2quat(euler(randi(length(uniqueGrains)),2:end)); %random orientation
fprintf(fileID,"        O: [%1.10f, %1.10f,%1.10f, %1.10f]\n",...
    quat(1),quat(2),quat(3),quat(4));

%print the rest of the properties
for i=1:length(uniqueGrains)
fprintf(fileID,"  - homogenization: SX\n");
fprintf(fileID,"    constituents:\n");

if material==1
    fprintf(fileID,"      - phase: 316L\n");
elseif material==2
     fprintf(fileID,"      - phase: IN625\n");
else
     fprintf(fileID,"      - phase: Aluminum\n");
end

fprintf(fileID,"        v: 1.0\n");
%fprintf(fileID,"        xi_0_sl: [%3.5e]\n",G*b*(HP_mult/avgDist+sqrt(disDensity))+PN_stress); % strengthening parameter

quat = eul2quat(euler(i,2:end));
    
fprintf(fileID,"        O: [%1.10f, %1.10f,%1.10f, %1.10f]\n",...
    quat(1),quat(2),quat(3),quat(4));
end
fclose(fileID);

movefile("material.txt","./"+workingDirectory+"/material.yaml"); %make the film a yaml
%% Write the grainID.geom file

%Downselect the full size of the microstructure
gridDownSelected = zeros([finalGridResolution(1),finalGridResolution(1),finalGridResolution(1)]);
gridDownSelected1 = imresize(gridPorous,"nearest",'OutputSize',finalGridResolution);

for i = 1:size(gridDownSelected,1)
gridDownSelected(:,i,:) = imresize(reshape(gridDownSelected1(:,i,:),[size(gridDownSelected1,1),size(gridDownSelected1,3)])...
    ,"nearest",'OutputSize',finalGridResolution);
end

%relabel the grain ID's to start from 1
gridRelabeled = gridDownSelected;
uniqueRelabeled = unique(gridRelabeled);

for i = 1:length(uniqueRelabeled) %start labeled from 0
    gridRelabeled(gridDownSelected==uniqueRelabeled(i))= i-1;
end

fileID = fopen('grainIDs.txt','w');
fprintf(fileID,'5	header\n');
fprintf(fileID,'grid	a %i	b %i	c %i \n',finalGridResolution(1), finalGridResolution(1), finalGridResolution(1));
fprintf(fileID,'size	x 1.0	y 1.0	z 1.0 \norigin	x 0.0	y 0.0	z 0.0\nhomogenization	1\n<!skip>\n');

for z = 1:finalGridResolution(1)
    matToWrite = gridRelabeled(:,:,z)+1;
    dlmwrite('grainIDs.txt',matToWrite(:)','delimiter',' ','-append');
end

%Convert grain ID file to the correct postfix (.geom,.yaml)
movefile("grainIDs.txt","./"+workingDirectory+"/grainIDs.geom");

%% Copy the submission scripts
copyfile("submit_simulation_v1.sh","./"+workingDirectory+"/submit_simulation_v1.sh");

end