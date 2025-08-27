gridSize = 500;
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

%% Extract the grain size distributions

%Compare to simulated data
scanSpeed02 = [];
scanSpeed03 = [];
scanSpeed04 = [];
scanSpeed05 = [];

numfiles = 240;
filePathOrig ='./2023_04_13_ParametricStudy_Steel';

for trialPlot= 0:numfiles-1
    load(filePathOrig+"/trial_"+string(round(trialPlot,0))+"/mat_dat.mat");
    
    %Calculate the grain size
    allGrainSizes=[];
    u = unique(grid);
    for i = 1:length(u)
        allGrainSizes=[allGrainSizes;sum(sum(grid==u(i)))];
    end
    allGrainSizes = 2*sqrt(allGrainSizes*((300e-6/gridSize)*10^6)^2 / 3.1415);

    if scanSpeed ==0.2
        scanSpeed02 = [scanSpeed02;allGrainSizes];
    elseif scanSpeed ==0.3
        scanSpeed03 = [scanSpeed03;allGrainSizes];
    elseif scanSpeed ==0.4
        scanSpeed04 = [scanSpeed04;allGrainSizes];
    elseif scanSpeed ==0.5
        scanSpeed05 = [scanSpeed05;allGrainSizes];
    end
end

bins = 1:0.5:60;

subplot(4,1,1)
histogram(scanSpeed02,bins,'Normalization', 'pdf');
xlim([0,50])

subplot(4,1,2)
histogram(scanSpeed03,bins,'Normalization', 'pdf');
xlim([0,50])

subplot(4,1,3)
histogram(scanSpeed04,bins,'Normalization', 'pdf');
xlim([0,50])

subplot(4,1,4)
histogram(scanSpeed05,bins,'Normalization', 'pdf');
xlim([0,50])

