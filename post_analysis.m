%% Post Analysis of a parametric study file
numFiles=240;

materialIndicator_Data = [];
scanSpeed_Data = [];
laserPower_Data = [];
crossHatch_Data = [];
heterogeneousShapeFactor_Data = [];
uniqueGrains = [];
asymmetryRatio = [];
grainSize=[];


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

    %Calculate the asymmetry ratio
    pAxes = findprincipleAxes(grid,xgrid,ygrid,[],cmap);
    pMax = max(pAxes(:,2:3),[],2);
    pMin = min(pAxes(:,2:3),[],2);

    asymmetryRatio=[asymmetryRatio;median(pMax./pMin)];

    %subplot(9,3,trialPlot+1)
    %imshow(grid+2,cmap)
    %title(sprintf("v=%.1f,P=%i,crossHatch=%.2f",[scanSpeed,laserPower,crossHatch]))
    %axis off

end

%grainSize = 300*300./uniqueGrains;
%grainSize = sqrt(grainSize/3.1415);

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
%grain asymmetry effects

figure
subplot(3,2,1)

gscatter(laserPower_Data,asymmetryRatio,crossHatch_Data)
ylabel("Grain Asymmetry Ratio [b/a]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,2)

gscatter(laserPower_Data,asymmetryRatio,scanSpeed_Data)
ylabel("Grain Asymmetry Ratio [b/a]")
xlabel("Laser Power [W]")
l=legend;
l.Title.String ="Scan Speed [m/s]";

subplot(3,2,3)

gscatter(scanSpeed_Data,asymmetryRatio,laserPower_Data)
ylabel("Grain Asymmetry Ratio [b/a]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,4)

gscatter(scanSpeed_Data,asymmetryRatio,crossHatch_Data)
ylabel("Grain Asymmetry Ratio [b/a]")
xlabel("Scan Speed [m/s]")
l=legend;
l.Title.String ="Cross Hatch";

subplot(3,2,5)

gscatter(crossHatch_Data,asymmetryRatio,laserPower_Data)
ylabel("Grain Asymmetry Ratio [b/a]")
xlabel("Cross Hatch")
l=legend;
l.Title.String ="Laser Power [W]";

subplot(3,2,6)

gscatter(crossHatch_Data,asymmetryRatio,scanSpeed_Data)
ylabel("Grain Asymmetry Ratio [b/a]")
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
t.asymmetry_ratio = asymmetryRatio;

tbstats1 = grpstats(t,["ScanSpeed","crossHatch"],["mean","std"],"DataVars",["grainSize"]);
tbstats2 = grpstats(t,["ScanSpeed","laserPower"],["mean","std"],"DataVars",["grainSize"]);
tbstats3 = grpstats(t,["crossHatch","laserPower"],["mean","std"],"DataVars",["grainSize"]);

tbstats4 = grpstats(t,["ScanSpeed","crossHatch"],["mean","std"],"DataVars",["asymmetry_ratio"]);
tbstats5 = grpstats(t,["ScanSpeed","laserPower"],["mean","std"],"DataVars",["asymmetry_ratio"]);
tbstats6 = grpstats(t,["crossHatch","laserPower"],["mean","std"],"DataVars",["asymmetry_ratio"]); 


%% Calculation of volumetric Energy density 

%Tm = 1084+273; %melting point for copper, K
Tm = 1375+273; %melting point for steel, K

%Pre-calculate the pool diameter based on the scan speed and laser power
uniqueScan = unique(scanSpeed_Data);
uniqueLaserPower = unique(laserPower_Data);

beamDiameter = zeros(length(uniqueScan)*length(uniqueLaserPower),3);
VED = []; %volumetric energy density

iterVal = 1;
for indScan = 1:length(uniqueScan)
    for indLaserPower = 1:length(uniqueLaserPower)
        fprintf("iteration %i\n",iterVal);
        tempData = heat_fluxFunc(uniqueScan(indScan)*1000,uniqueLaserPower(indLaserPower),Tm);
        beamD = tempData.poolDiameter; %in m
        beamDiameter(iterVal,1) = uniqueScan(indScan);
        beamDiameter(iterVal,2) = uniqueLaserPower(indLaserPower);
        beamDiameter(iterVal,3) = beamD;
        iterVal=iterVal+1;
    end
end

%Calculate the VED for each set of processing parameters
for i = 1:height(t)
    beamD = beamDiameter(beamDiameter(:,1)==scanSpeed_Data(i) & beamDiameter(:,2)==laserPower_Data(i) , 3);
    VED_val = laserPower_Data(i)/(2*scanSpeed_Data(i)*beamD^2*crossHatch_Data(i));
    VED = [VED;VED_val];
end


%% Linear Fitting - asymmetry

%linear fitting - grain size 
grainAssymetry = @(k, data) k(1) + k(2).*data(:,1)+ k(3).*data(:,2) + k(4).*data(:,3);
k0 = rand(1,4);
xdata = [t.laserPower(:)/max(t.laserPower(:)),...
    t.ScanSpeed(:)/max(t.ScanSpeed(:)),...
    t.crossHatch(:)/max(t.crossHatch(:))];
ydata = (t.asymmetry_ratio-min(t.asymmetry_ratio))/(max(t.asymmetry_ratio)-min(t.asymmetry_ratio));
[B_grainAssy_lin,R,J,CovB]  = nlinfit(xdata, ydata, grainAssymetry,k0); 
