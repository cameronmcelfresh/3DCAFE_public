totalSims = 3;

subplot(3,totalSims,1);

for plotCounter = 1:totalSims
    
fileName = '/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_IN625_beamAdjusted_11_case_' ...
    + string(plotCounter) + '/workspace_' + string(plotCounter) + '.mat';
load(fileName);

% grainDetails_zNorm = grainSizeAsymmetry(grid(:,:,round(params{paramIter,7}*90/100)), dx, xgrid,ygrid);
% grainDetails_meltCenter = grainSizeAsymmetry(grid(:,round(params{paramIter,7}*52/100),round(20*params{paramIter,7}/100):end), dx, xgrid,ygrid);
% grainDetails_meltEdge = grainSizeAsymmetry(grid(:,round(63*params{paramIter,7}/100),round(20*params{paramIter,7}/100):end), dx, xgrid,ygrid);

grainDetails_zNorm = grainSizeAsymmetry(grid(:,:,90), dx, xgrid,ygrid);
grainDetails_meltCenter = grainSizeAsymmetry(grid(:,52,20:end), dx, xgrid,ygrid);
grainDetails_meltEdge = grainSizeAsymmetry(grid(:,63,20:end), dx, xgrid,ygrid);

subplot(3,totalSims,plotCounter);
histogram(grainDetails_zNorm(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(0.25,0.25,sprintf("mean = %.2f",mean(grainDetails_zNorm(:,2)*10^6)));
title("zNorm")

subplot(3,totalSims,totalSims+plotCounter);
histogram(grainDetails_meltCenter(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(0.25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltCenter(:,2)*10^6)));
title("meltCenter")

subplot(3,totalSims,2*totalSims+plotCounter);
histogram(grainDetails_meltEdge(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(0.25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltEdge(:,2)*10^6)));
title("meltEdge")

% %Asymmetry
% subplot(3,totalSims,plotCounter);
% histogram(grainDetails_zNorm(:,3),(1:0.5:6),'Normalization','probability');
% ylim([0,0.4])
% text(0.25,0.25,sprintf("mean = %.2f",mean(grainDetails_zNorm(:,3))));
% title("zNorm")
% 
% subplot(3,totalSims,totalSims+plotCounter);
% histogram(grainDetails_meltCenter(:,3),(1:0.5:6),'Normalization','probability');
% ylim([0,0.4])
% text(0.25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltCenter(:,3))));
% title("meltCenter")
% 
% subplot(3,totalSims,2*totalSims+plotCounter);
% histogram(grainDetails_meltEdge(:,3),(1:0.5:6),'Normalization','probability');
% ylim([0,0.4])
% text(0.25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltEdge(:,3))));
% title("meltEdge")


end


