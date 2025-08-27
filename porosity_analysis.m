
% fileDirs = {{'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v2_2/workspace_2.mat'},...
%     {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v2_1/workspace_1.mat'},...
%     {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_3/workspace_3.mat'},...
%     {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_2/workspace_2.mat'},...
%     {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_1/workspace_1.mat'},...
%     {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_5/workspace_5.mat'},...
%     {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_4/workspace_4.mat'}};

fileDirs = {{'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v3_1/workspace_1.mat'},...
    {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v2_1/workspace_1.mat'},...
    {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_3/workspace_3.mat'},...
    {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_2/workspace_2.mat'},...
    {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_1/workspace_1.mat'},...
    {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_5/workspace_5.mat'},...
    {'/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/temp_case_v1_4/workspace_4.mat'}};


porosity = zeros([1,7]);
stdError = zeros([1,7]);
VED = [26,32,39,46,65,88,117];

for iter_val = 1:7
    %load(fileDir + sprintf("temp_case_v1_%i/workspace_%i.mat",iter_val,iter_val));
    file=fileDirs{iter_val};
    load(file{1});

    porositySamples =[];
    for j = 1:5
        gridPorous = createPorosity(grid,xgrid,ygrid,zgrid,Tm*1.25,maxTemps,25);
        porositySamples = [porositySamples;sum(sum(sum(gridPorous==0)))/numel(gridPorous)*100;];
    end

    porosity(iter_val) = mean(porositySamples);
    stdError(iter_val) = std(porositySamples);
end

%plot(porosity)

errorbar(porosity,stdError)
figure
errorbar(VED,porosity,stdError)
ylabel("Porosity")
xlabel("VED")