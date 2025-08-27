
load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_316L_5_case_1/workspace_1.mat')
grainDetails_zNorm1 = grainSizeAsymmetry(grid(:,:,60), dx, xgrid,ygrid);
grainDetails_meltCenter1 = grainSizeAsymmetry(grid(:,40,:), dx, xgrid,ygrid);
grainDetails_meltEdge1 = grainSizeAsymmetry(grid(:,63,:), dx, xgrid,ygrid);

load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_316L_5_case_2/workspace_2.mat')
grainDetails_zNorm2 = grainSizeAsymmetry(grid(:,:,60), dx, xgrid,ygrid);
grainDetails_meltCenter2 = grainSizeAsymmetry(grid(:,40,:), dx, xgrid,ygrid);
grainDetails_meltEdge2 = grainSizeAsymmetry(grid(:,63,:), dx, xgrid,ygrid);

load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_316L_5_case_3/workspace_3.mat')
grainDetails_zNorm3 = grainSizeAsymmetry(grid(:,:,70), dx, xgrid,ygrid);
grainDetails_meltCenter3 = grainSizeAsymmetry(grid(:,40,:), dx, xgrid,ygrid);
grainDetails_meltEdge3 = grainSizeAsymmetry(grid(:,63,:), dx, xgrid,ygrid);


subplot(3,3,1);
histogram(grainDetails_zNorm1(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_zNorm1(:,2)*10^6)));
title("zNorm")
subplot(3,3,2);
histogram(grainDetails_zNorm2(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_zNorm2(:,2)*10^6)));
title("zNorm")
subplot(3,3,3);
histogram(grainDetails_zNorm3(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_zNorm3(:,2)*10^6)));
title("zNorm")

subplot(3,3,4);
histogram(grainDetails_meltCenter1(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltCenter1(:,2)*10^6)));
title("meltCenter")
subplot(3,3,5);
histogram(grainDetails_meltCenter2(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltCenter2(:,2)*10^6)));
title("meltCenter")
subplot(3,3,6);
histogram(grainDetails_meltCenter3(:,2)*10^6,(5:5:50),'Normalization','probability');
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltCenter3(:,2)*10^6)));
ylim([0,0.4])
title("meltCenter")


subplot(3,3,7);
histogram(grainDetails_meltEdge1(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltEdge1(:,2)*10^6)));
title("meltEdge")
subplot(3,3,8);
histogram(grainDetails_meltEdge2(:,2)*10^6,(5:5:50),'Normalization','probability');
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltEdge2(:,2)*10^6)));
ylim([0,0.4])
title("meltEdge")
subplot(3,3,9);
histogram(grainDetails_meltEdge3(:,2)*10^6,(5:5:50),'Normalization','probability');
text(25,0.25,sprintf("mean = %.2f",mean(grainDetails_meltEdge3(:,2)*10^6)));
ylim([0,0.4])
title("meltEdge")


