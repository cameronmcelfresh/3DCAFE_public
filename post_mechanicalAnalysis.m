%% Post Analysis of a parametric study file

numFiles=240;

materialIndicator_Data = [];
scanSpeed_Data = [];
laserPower_Data = [];
crossHatch_Data = [];
heterogeneousShapeFactor_Data = [];
uniqueGrains = [];

yieldStrength = [];
grainSize=[];
hardeningRate = [];

%filePathOrig="./2023_03_15_ParametricStudy";
%filePathOrig='./2023_03_18_ParametricStudy_Cu';
filePathOrig="./2023_04_07_ParametricStudy_Steel";
%filePathOrig="./2023_04_13_ParametricStudy_Steel";

for trialPlot= 0:numFiles-1

    load(filePathOrig+"/trial_"+string(round(trialPlot,0))+"/mat_dat.mat");

    materialIndicator_Data = [materialIndicator_Data;materialIndicator];
    scanSpeed_Data = [scanSpeed_Data;scanSpeed];
    laserPower_Data = [laserPower_Data;laserPower];
    crossHatch_Data = [crossHatch_Data;crossHatch];
    heterogeneousShapeFactor_Data = [heterogeneousShapeFactor_Data;het_shapeFactor];
    uniqueGrains = [uniqueGrains;length(unique(grid))];

    %Calculate the grain size
    allGrainSizes=[];
    u = unique(grid);
    for i = 1:length(u)
        allGrainSizes=[allGrainSizes;sum(sum(grid==u(i)))];
    end
    allGrainSizes = sqrt(allGrainSizes*((realGridSize/gridSize)*10^6)^2 / 3.1415);
    grainSize = [grainSize;mean(allGrainSizes)*2];

    %Calculate the yield strength
    [yieldStrength_val,hardeningRate_val] = calc_strengthData(filePathOrig+"/trial_"+string(round(trialPlot,0))+"/stress_strain.csv");

    yieldStrength=[yieldStrength;yieldStrength_val];
    hardeningRate = [hardeningRate;hardeningRate_val];
end

%% Plot results

%grain size effects
figure
subplot(3,2,1)

gscatter(laserPower_Data,grainSize,crossHatch_Data)
ylabel("Grain Size [um]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,2)

gscatter(laserPower_Data,grainSize,scanSpeed_Data)
ylabel("Grain Size [um]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Scan Speed [m/s]";

subplot(3,2,3)

gscatter(scanSpeed_Data,grainSize,laserPower_Data)
ylabel("Grain Size [um]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,4)

gscatter(scanSpeed_Data,grainSize,crossHatch_Data)
ylabel("Grain Size [um]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,5)

gscatter(crossHatch_Data,grainSize,laserPower_Data)
ylabel("Grain Size [um]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,6)




gscatter(crossHatch_Data,grainSize,scanSpeed_Data)
ylabel("Grain Size [um]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Scan Speed [m/s]";

%%%
%Yield Strength effects

figure
subplot(3,2,1)

gscatter(laserPower_Data,yieldStrength,crossHatch_Data)
ylabel("Yield Strength [MPa]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,2)

gscatter(laserPower_Data,yieldStrength,scanSpeed_Data)
ylabel("Yield Strength [MPa]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Scan Speed [m/s]";

subplot(3,2,3)

gscatter(scanSpeed_Data,yieldStrength,laserPower_Data)
ylabel("Yield Strength [MPa]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,4)

gscatter(scanSpeed_Data,yieldStrength,crossHatch_Data)
ylabel("Yield Strength [MPa]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,5)

gscatter(crossHatch_Data,yieldStrength,laserPower_Data)
ylabel("Yield Strength [MPa]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,6)

gscatter(crossHatch_Data,yieldStrength,scanSpeed_Data)
ylabel("Yield Strength [MPa]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Scan Speed [m/s]";

%%%
%Hardening Rate Effects

figure
subplot(3,2,1)

gscatter(laserPower_Data,hardeningRate,crossHatch_Data)
ylabel("Hardening Rate [MPa/m/m]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,2)

gscatter(laserPower_Data,hardeningRate,scanSpeed_Data)
ylabel("Hardening Rate [MPa/m/m]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Scan Speed [m/s]";

subplot(3,2,3)

gscatter(scanSpeed_Data,hardeningRate,laserPower_Data)
ylabel("Hardening Rate [MPa/m/m]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,4)

gscatter(scanSpeed_Data,hardeningRate,crossHatch_Data)
ylabel("Hardening Rate [MPa/m/m]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,5)

gscatter(crossHatch_Data,hardeningRate,laserPower_Data)
ylabel("Hardening Rate [MPa/m/m]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,6)

gscatter(crossHatch_Data,hardeningRate,scanSpeed_Data)
ylabel("Hardening Rate [MPa/m/m]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Scan Speed [m/s]";

%% Make output data tables

t = table();
t.ScanSpeed = scanSpeed_Data;
t.laserPower = laserPower_Data;
t.crossHatch = crossHatch_Data;
t.uniqueGrains = uniqueGrains;
t.grainSize = grainSize;
t.yieldStrength = yieldStrength;
t.hardeningRate = hardeningRate;

tbstats1 = grpstats(t,["ScanSpeed","crossHatch"],["mean","std"],"DataVars",["grainSize"]);
tbstats2 = grpstats(t,["ScanSpeed","laserPower"],["mean","std"],"DataVars",["grainSize"]);
tbstats3 = grpstats(t,["crossHatch","laserPower"],["mean","std"],"DataVars",["grainSize"]);

tbstats4 = grpstats(t,["ScanSpeed","crossHatch"],["mean","std"],"DataVars",["yieldStrength"]);
tbstats5 = grpstats(t,["ScanSpeed","laserPower"],["mean","std"],"DataVars",["yieldStrength"]);
tbstats6 = grpstats(t,["crossHatch","laserPower"],["mean","std"],"DataVars",["yieldStrength"]);

tbstats7 = grpstats(t,["ScanSpeed","crossHatch"],["mean","std"],"DataVars",["hardeningRate"]);
tbstats8 = grpstats(t,["ScanSpeed","laserPower"],["mean","std"],"DataVars",["hardeningRate"]);
tbstats9 = grpstats(t,["crossHatch","laserPower"],["mean","std"],"DataVars",["hardeningRate"]); 


%% Non-linear Variable Fitting

%non-linear fitting, yield strength
mechF = @(k, data) k(1).*data(:,1).^k(2) + k(3).*data(:,2).^k(4) + k(5).*data(:,3).^k(6);
k0 = rand(1,6);
xdata = [t.laserPower(:)/max(t.laserPower(:)),...
    t.ScanSpeed(:)/max(t.ScanSpeed(:)),...
    t.crossHatch(:)/max(t.crossHatch(:))];
ydata = t.yieldStrength;
[B_yield,R,J,CovB]  = nlinfit(xdata, ydata, mechF,k0); % <-- B_yield has all of the constants
%BCI_yield = nlparci(B_yield,R,'covar',CovB);

ydata = t.hardeningRate;
[B_hardening,R,J,CovB]  = nlinfit(xdata, ydata, mechF,k0); % <-- B_hardening has all of the constants
%BCI_hardening = nlparci(B_yield,R,'covar',CovB);


%linear fitting
mechF = @(k, data) k(1) + k(2).*data(:,1)+ k(3).*data(:,2) + k(4).*data(:,3);
k0 = rand(1,4);
xdata = [t.laserPower(:)/max(t.laserPower(:)),...
    t.ScanSpeed(:)/max(t.ScanSpeed(:)),...
    t.crossHatch(:)/max(t.crossHatch(:))];
ydata = (t.yieldStrength-min(t.yieldStrength))/(max(t.yieldStrength)-min(t.yieldStrength));
[B_yield_lin,R,J,CovB]  = nlinfit(xdata, ydata, mechF,k0); % <-- B_yield_lin has all of the constants
%BCI_yield = nlparci(B_yield,R,'covar',CovB);

ydata = (t.hardeningRate-min(t.hardeningRate))/(max(t.hardeningRate)-min(t.hardeningRate));
[B_hardening_lin,R,J,CovB]  = nlinfit(xdata, ydata, mechF,k0); % <-- B_hardening_lin has all of the constants
%BCI_hardening = nlparci(B_yield,R,'covar',CovB);

%linear fitting - grain size 
grainF = @(k, data) k(1) + k(2).*data(:,1)+ k(3).*data(:,2) + k(4).*data(:,3);
k0 = rand(1,4);
xdata = [t.laserPower(:)/max(t.laserPower(:)),...
    t.ScanSpeed(:)/max(t.ScanSpeed(:)),...
    t.crossHatch(:)/max(t.crossHatch(:))];
ydata = (t.grainSize-min(t.grainSize))/(max(t.grainSize)-min(t.grainSize));
[B_grainSize_lin,R,J,CovB]  = nlinfit(xdata, ydata, grainF,k0); % <-- B_yield_lin has all of the constants

