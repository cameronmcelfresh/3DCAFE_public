gridSize = 800;
[xgrid,ygrid] = meshgrid((1:gridSize),(1:gridSize));

numStructures= 10;
asymRatio = [];

for i = 1:numStructures
    [grid,~] = growEquiaxed(gridSize,500);
    pAxes = findprincipleAxes(grid,xgrid,ygrid,[],[]);
    asymRatio = [asymRatio;max(pAxes(:,2:3),[],2)./min(pAxes(:,2:3),[],2)];
end

h=histogram(asymRatio,'Normalization', 'pdf');
b = h.BinEdges;
counts = h.BinCounts;

%% Extract the grain asymmety ratio
targetScanSpeed = 0.2;

%Compare to simulated data
asymmetry100 = [];
asymmetry150 = [];
asymmetry200 = [];
asymmetry250= [];
asymmetry300 = [];

numfiles = 240;
filePathOrig ='./2023_03_18_ParametricStudy_Cu';

for trialPlot= 0:numfiles-1
    load(filePathOrig+"/trial_"+string(round(trialPlot,0))+"/mat_dat.mat");
    
    if scanSpeed ~= targetScanSpeed
        continue;
    end
    %Calculate the asymmetry ratio
    pAxes = findprincipleAxes(grid,xgrid,ygrid,[],cmap);
    pMax = max(pAxes(:,2:3),[],2);
    pMin = min(pAxes(:,2:3),[],2);

    if laserPower ==100
        asymmetry100= [asymmetry100;max(pAxes(:,2:3),[],2)./min(pAxes(:,2:3),[],2)];
    elseif laserPower ==150
        asymmetry150= [asymmetry150;max(pAxes(:,2:3),[],2)./min(pAxes(:,2:3),[],2)];        
    elseif laserPower ==200
        asymmetry200= [asymmetry200;max(pAxes(:,2:3),[],2)./min(pAxes(:,2:3),[],2)];
    elseif laserPower ==250
        asymmetry250= [asymmetry250;max(pAxes(:,2:3),[],2)./min(pAxes(:,2:3),[],2)];
    elseif laserPower ==300
        asymmetry300= [asymmetry300;max(pAxes(:,2:3),[],2)./min(pAxes(:,2:3),[],2)];
    end
end

asymmetry100(asymmetry100>1e3)=[];
asymmetry150(asymmetry150>1e3)=[];
asymmetry200(asymmetry200>1e3)=[];
asymmetry250(asymmetry250>1e3)=[];
asymmetry300(asymmetry300>1e3)=[];

bins = 1:0.05:4;

subplot(6,1,1)
histogram(asymRatio,bins,'Normalization', 'pdf');
xlim([1,3])
ylim([0,1.5])

subplot(6,1,2)
histogram(asymmetry100,bins,'Normalization', 'pdf');
xlim([1,3])
ylim([0,1.5])

subplot(6,1,3)
histogram(asymmetry150,bins,'Normalization', 'pdf');
xlim([1,3])
ylim([0,1.5])

subplot(6,1,4)
histogram(asymmetry200,bins,'Normalization', 'pdf');
xlim([1,3])
ylim([0,1.5])

subplot(6,1,5)
histogram(asymmetry250,bins,'Normalization', 'pdf');
xlim([1,3])
ylim([0,1.5])

subplot(6,1,6)
histogram(asymmetry300,bins,'Normalization', 'pdf');
xlim([1,3])
ylim([0,1.5])

legend
