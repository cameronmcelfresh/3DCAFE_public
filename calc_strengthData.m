function [yieldStrength,hardeningRate_avg] = calc_strengthData(filePath)
%calc_strengthData Function to calculate the yield strength and hardening
%rate from a CPFE simulation

%Read in the data
mat = readmatrix(filePath);

%% Find the yield strength with 0.002 offset

%Oversample with interpolation
strainOverSample = linspace(0,max(mat(:,1)),100);
stressOverSample = interp1(mat(:,1),mat(:,2),strainOverSample);

%Manually calculate the slope (modulus)
slope = (mat(2,2)-mat(1,2)) / (mat(2,1)-mat(1,1));

%Find where the two curves intersect
[strainIntersect,stressIntersect] = polyxpoly(strainOverSample,stressOverSample,...
    strainOverSample+0.002,strainOverSample*slope);

% %Plot the 0.002 offset
% figure
% plot(strainOverSample,stressOverSample/10^6);
% hold on
% plot(strainOverSample+0.002,(strainOverSample*slope)/10^6);
% 
% % %Plot the intersection point
% scatter(strainIntersect,stressIntersect/10^6,20,'r','filled');
% xlabel("Strain");
% ylabel("Stress [MPa]");
% ylim([0,max(stressOverSample)/10^6*1.05]);

yieldStrength = stressIntersect/10^6;

%% Calculate the hardening Rate

%Find the plastic portion of the stress and strain
plasticStress = mat(mat(:,2)>yieldStrength,2)/10^6; % in MPa
plasticStrain = mat(mat(:,2)>yieldStrength,1);

%Calcualte the 1 dimensional gradient for strain hardening
Fstress = gradient(plasticStress);
Fstrain = gradient(plasticStrain);

hardeningRate = Fstress./Fstrain; %calculate the hardening rate

% %Plot
% figure
% plot(plasticStrain,hardeningRate);
% xlabel("Plastic Strain");
% ylabel("Hardening Rate");

hardeningRate_avg=median(hardeningRate);
end
