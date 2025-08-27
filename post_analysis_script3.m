
load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_rotate_1/workspace_1.mat')
grainDetails_zNorm1 = grainSizeAsymmetry(grid(:,:,90), dx, xgrid,ygrid);
grainDetails_meltCenter1 = grainSizeAsymmetry(grid(:,50,:), dx, xgrid,ygrid);
grainDetails_meltEdge1 = grainSizeAsymmetry(grid(:,45,:), dx, xgrid,ygrid);

load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_rotate_2/workspace_2.mat')
grainDetails_zNorm2 = grainSizeAsymmetry(grid(:,:,90), dx, xgrid,ygrid);
grainDetails_meltCenter2 = grainSizeAsymmetry(grid(:,50,:), dx, xgrid,ygrid);
grainDetails_meltEdge2 = grainSizeAsymmetry(grid(:,45,:), dx, xgrid,ygrid);


subplot(3,2,1);
histogram(grainDetails_zNorm1(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("median = %.2f",mean(grainDetails_zNorm1(:,2)*10^6)));
title("zNorm")
subplot(3,2,2);
histogram(grainDetails_zNorm2(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("median = %.2f",mean(grainDetails_zNorm2(:,2)*10^6)));
title("zNorm")

subplot(3,2,3);
histogram(grainDetails_meltCenter1(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("median = %.2f",mean(grainDetails_meltCenter1(:,2)*10^6)));
title("meltCenter")
subplot(3,2,4);
histogram(grainDetails_meltCenter2(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("median = %.2f",mean(grainDetails_meltCenter2(:,2)*10^6)));
title("meltCenter")


subplot(3,2,5);
histogram(grainDetails_meltEdge1(:,2)*10^6,(5:5:50),'Normalization','probability');
ylim([0,0.4])
text(25,0.25,sprintf("median = %.2f",mean(grainDetails_meltEdge1(:,2)*10^6)));
title("meltEdge")
subplot(3,2,6);
histogram(grainDetails_meltEdge2(:,2)*10^6,(5:5:50),'Normalization','probability');
text(25,0.25,sprintf("median = %.2f",mean(grainDetails_meltEdge2(:,2)*10^6)));
ylim([0,0.4])
title("meltEdge")


%% write the text files
% titles = {'grainSize_um','prob'};
% 
% figure
% h=histogram(grainDetails_zNorm1(:,2)*10^6,(5:5:50),'Normalization','probability');
% case1_zNorm = [h.Values,0];
% bins=[h.BinEdges];
% h=histogram(grainDetails_zNorm2(:,2)*10^6,(5:5:50),'Normalization','probability');
% case2_zNorm = [h.Values,0];
% h=histogram(grainDetails_zNorm3(:,2)*10^6,(5:5:50),'Normalization','probability');
% case3_zNorm = [h.Values,0];
% h=histogram(grainDetails_zNorm4(:,2)*10^6,(5:5:50),'Normalization','probability');
% case4_zNorm = [h.Values,0];
% 
% h=histogram(grainDetails_meltCenter1(:,2)*10^6,(5:5:50),'Normalization','probability');
% case1_meltCenter = [h.Values,0];
% h=histogram(grainDetails_meltCenter2(:,2)*10^6,(5:5:50),'Normalization','probability');
% case2_meltCenter = [h.Values,0];
% h=histogram(grainDetails_meltCenter3(:,2)*10^6,(5:5:50),'Normalization','probability');
% case3_meltCenter = [h.Values,0];
% h=histogram(grainDetails_meltCenter4(:,2)*10^6,(5:5:50),'Normalization','probability');
% case4_meltCenter = [h.Values,0];
% 
% h=histogram(grainDetails_meltEdge1(:,2)*10^6,(5:5:50),'Normalization','probability');
% case1_meltEdge = [h.Values,0];
% h=histogram(grainDetails_meltEdge2(:,2)*10^6,(5:5:50),'Normalization','probability');
% case2_meltEdge = [h.Values,0];
% h=histogram(grainDetails_meltEdge3(:,2)*10^6,(5:5:50),'Normalization','probability');
% case3_meltEdge = [h.Values,0];
% h=histogram(grainDetails_meltEdge4(:,2)*10^6,(5:5:50),'Normalization','probability');
% case4_meltEdge = [h.Values,0];
% 
% C = [titles;num2cell([bins(:),case1_zNorm(:)])];
% writecell(C,'case1_zNorm_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case2_zNorm(:)])];
% writecell(C,'case2_zNorm_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case3_zNorm(:)])];
% writecell(C,'case3_zNorm_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case4_zNorm(:)])];
% writecell(C,'case4_zNorm_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case1_meltCenter(:)])];
% writecell(C,'case1_meltCenter_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case2_meltCenter(:)])];
% writecell(C,'case2_meltCenter_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case3_meltCenter(:)])];
% writecell(C,'case3_meltCenter_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case4_meltCenter(:)])];
% writecell(C,'case4_meltCenter_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case1_meltEdge(:)])];
% writecell(C,'case1_meltEdge_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case2_meltEdge(:)])];
% writecell(C,'case2_meltEdge_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case3_meltEdge(:)])];
% writecell(C,'case3_meltEdge_GS.txt','Delimiter','space')
% 
% C = [titles;num2cell([bins(:),case4_meltEdge(:)])];
% writecell(C,'case4_meltEdge_GS.txt','Delimiter','space')

%% aspect ratio
figure

subplot(3,2,1);
histogram(grainDetails_zNorm1(:,3),(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",median(grainDetails_zNorm1(:,3))));
ylim([0,0.5])
title("zNorm")
subplot(3,2,2);
histogram(grainDetails_zNorm2(:,3),(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_zNorm2(:,3))));
ylim([0,0.5])
title("zNorm")

subplot(3,2,3);
histogram(grainDetails_meltCenter1(:,3),(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltCenter1(:,3))));
ylim([0,0.5])
title("meltCenter")
subplot(3,2,4);
histogram(grainDetails_meltCenter2(:,3),(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltCenter2(:,3))));
ylim([0,0.5])
title("meltCenter")

subplot(3,2,5);
histogram(grainDetails_meltEdge1(:,3),(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltEdge1(:,3))));
ylim([0,0.5])
title("meltEdge")
subplot(3,2,6);
histogram(grainDetails_meltEdge2(:,3),(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltEdge2(:,3))));
ylim([0,0.5])
title("meltEdge")

%volume fraction
figure

subplot(3,2,1);
histogram(grainDetails_zNorm1(:,3).^2,(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",median(grainDetails_zNorm1(:,3))));
ylim([0,0.5])
title("zNorm - area fraction")
subplot(3,2,2);
histogram(grainDetails_zNorm2(:,3).^2,(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_zNorm2(:,3))));
ylim([0,0.5])
title("zNorm - area fraction")


subplot(3,2,3);
histogram(grainDetails_meltCenter1(:,3).^2,(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltCenter1(:,3))));
ylim([0,0.5])
title("meltCenter - area fraction")
subplot(3,2,4);
histogram(grainDetails_meltCenter2(:,3).^2,(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltCenter2(:,3))));
ylim([0,0.5])
title("meltCenter - area fraction")


subplot(3,2,5);
histogram(grainDetails_meltEdge1(:,3).^2,(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltEdge1(:,3))));
ylim([0,0.5])
title("meltEdge - area fraction")
subplot(3,2,6);
histogram(grainDetails_meltEdge2(:,3).^2,(1:0.5:6),'Normalization','probability');
text(3,0.25,sprintf("median = %.2f",mean(grainDetails_meltEdge2(:,3))));
ylim([0,0.5])
title("meltEdge - area fraction")

%%
%% Texture Analysis
figure
mtexFig = newMtexFigure('layout',[1,2]);

cs = crystalSymmetry('cubic');

nextAxis(1,1)

load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_rotate_1/workspace_1.mat')
gridEuler1 = grid(xgrid>50 & xgrid<60 & zgrid>30);
g = unique(gridEuler1(:));
g=g(g>0);
mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF - 1")
setColorRange([0.2,3])
t = norm(odf)^2

nextAxis(1,2)

load('/Users/cameronmcelfresh/Documents/UCLA Research/Additive - 3D CPFE/MATLAB_CODE_3D/sim_rotate_2/workspace_2.mat')
gridEuler1 = grid(xgrid>50 & xgrid<60 & zgrid>30);
g = unique(gridEuler1(:));
g=g(g>0);
mod1 = rotation.byEuler(euler(g-1000,2),euler(g-1000,3),euler(g-1000,4),cs);
odf = unimodalODF(mod1);
plotIPDF(odf,vector3d(0,0,1),'antipodal')
title("001 - IPF, 2")
setColorRange([0.2,3])
t = norm(odf)^2

%
% mtexFig = newMtexFigure('layout',[2,3]);
% 
% plot(xvector,'upper')
% 
% nextAxis
% 
% plot(zvector,'upper')
% 
% nextAxis(2,3)