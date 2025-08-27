function [C,R,y_normPDF,Rs] = spherePacking(boxSize,R_mean,R_STD,maxIters)
%spherePacking Packs Spheres into a 3D box
%   Detailed explanation goes here

C = [];
R = [];

queryRadius = R_mean-10*R_STD:R_mean/100:R_mean+10*R_STD;
y_normPDF= normpdf(queryRadius, R_mean,R_STD);
y_normPDF=y_normPDF/max(y_normPDF);
Rs=queryRadius;
for i=1:maxIters
    newPoint = rand([1,3])*boxSize;

    if i==1
        C=[newPoint];
        R = [R_mean];
        continue;
    end

    %dists = C-newPoint;
    %dists = sqrt(C(1,:).^2+C(2,:).^2+C(:,3).^2);
    dists = pdist2(C,newPoint);
    [minDistance,minIndex] = min(dists);
    newRadius = minDistance-R(minIndex);
    probAccept = interp1(queryRadius,y_normPDF,newRadius); %sample the pdf to see the probability of accepting
    
    if probAccept>rand()
        C=[C;newPoint];
        R=[R;newRadius];
    end
end

% %plot
% figure
% for i = 1:length(R)
%     [x,y,z] = sphere;
%     radius = R(i);
%     x = x * radius;
%     y = y * radius;
%     z = z * radius;
%     % Translate sphere to new location.
%     % Plot as surface.
%     hh=surf(x+C(i,1),y+C(i,2),z+C(i,3));
%     set(hh,'EdgeColor','None','FaceAlpha',0.3,'FaceColor','b')
%     hold on
% end
% % Label axes.
% xlabel('X', 'FontSize', 20);
% ylabel('Y', 'FontSize', 20);
% zlabel('Z', 'FontSize', 20);
% 
%Estimate packing efficiency
in = 0;
totalPoints = maxIters;

for i = 1:maxIters
    newPoint = rand([1,3])*boxSize;
    dists = pdist2(C,newPoint);
    if any(dists<R)
        in=in+1;
    end
end

fprintf("\nEstimated volume infill: %.2f%%\n", in/totalPoints*100);

end