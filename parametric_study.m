%% Script to run a parametric study of 2D AM grain simulation

%Assign the file path
%filePath="./2023_04_29_ParametricStudy_Test";
filePath="./2023_05_01_paramTest";

%% Assign the variables to cycle through during the parametric study
% materialIndicator_array = [1];
% scanSpeed_array = [0.2,0.3,0.4,0.5];
% laserPower_array = [100,150,200,250,300];
% crossHatch_array = [0.1,0.2,0.3,0.4];
% heterogeneousShapeFactor_array = [0.5];
% duplicates = 3;

materialIndicator_array = [1];
scanSpeed_array = [0.3,0.4];
laserPower_array = [100];
crossHatch_array = [0.1];
heterogeneousShapeFactor_array = [0.5];
duplicates = 1;

restartIter = 0; %iteration to restart at

%% Run the parametric study
mainSolidification_CA_v1_function(filePath,...
    materialIndicator_array,... %array of indicators for the material type
    scanSpeed_array,... %array of scan speeds to run
    laserPower_array,... %array of laser powers to run
    crossHatch_array,... %array of the crosshatch spacing
    heterogeneousShapeFactor_array,... %array of shape factors
    duplicates,... %number of duplicates to run
    restartIter);  %iteration to restart at

%% Run post-processing analysis

