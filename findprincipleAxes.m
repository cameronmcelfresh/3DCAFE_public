function [axes] = findprincipleAxes(grid,xgrid,ygrid,dx,cmap)
%FINDPRINCIPLEAXES Function to find the principle axes for each unique
%grain

axes = [];
%column 1 = grain ID
%column 2 = x component of 1st principle axis
%column 3 = y component of 1st principle axis
%column 4 = x component of 2nd principle axis
%column 5 = y component of 2nd principle axis
%column 6 = rotation angle of tilted ellipse
%column 7 = latent of 1st principle axis
%column 8 = latent of 2nd principle axis
axes = zeros([length(unique(grid)),10]);

uniqueGrains = unique(grid); %find the unique grain identifiers

%For plotting
%imshow(grid,cmap)
%axis off
%hold on

for grainIndex = 1:length(uniqueGrains)
    %fprintf("Grain "+uniqueGrains(grainIndex)+"\n");
    
    %Find grain ID
    grainID = uniqueGrains(grainIndex);
    
    %Find first principle axis
    x = xgrid(grid==grainID);
    y = ygrid(grid==grainID);
    XY = [x(:),y(:)];
    mu = mean(XY);

    if length(XY)<8 %if less than 5 points are found then an ellipse cannot be fit
        axes(grainIndex,1)=uniqueGrains(grainIndex);
        axes(grainIndex,2)= sqrt(length(XY)/3.1415) ;
        axes(grainIndex,3)=sqrt(length(XY)/3.1415);
        axes(grainIndex,4)=0;
        axes(grainIndex,5)=0;
        axes(grainIndex,6)=0;
        axes(grainIndex,7)=0;
        axes(grainIndex,8)=0;        
        axes(grainIndex,9)=mu(1);
        
        if length(mu)==1
            axes(grainIndex,10)=mu(1); 
        else
            axes(grainIndex,10)=mu(2);      
        end
        continue;
    end

    %[princpAxes,~,latent] = pca(XY); %use pca to find the principle axes and latent
    %Find the center
    %PCA Attempt
    %princpAxes(1,:)=princpAxes(1,:)*latent(1)/100;
    %princpAxes(2,:)=princpAxes(2,:)*latent(2)/100;

    %Fit Ellipse
    ellipse_t = fit_ellipse(XY(:,1),XY(:,2) );

    if isempty(ellipse_t) %if the rotation matrix is empty then too few points were found for a given grain
        axes(grainIndex,1)=uniqueGrains(grainIndex);
        axes(grainIndex,2)= sqrt(length(XY)/3.1415) ;
        axes(grainIndex,3)=sqrt(length(XY)/3.1415);
        axes(grainIndex,4)=0;
        axes(grainIndex,5)=0;
        axes(grainIndex,6)=0;
        axes(grainIndex,7)=0;
        axes(grainIndex,8)=0;        
        axes(grainIndex,9)=mu(1);
        axes(grainIndex,10)=mu(2);          
        continue;
    end

    R = [cos(ellipse_t.phi), -sin(ellipse_t.phi);
        sin(ellipse_t.phi), cos(ellipse_t.phi)];
    
    if isempty(R) %if the rotation matrix is empty then too few points were found for a given grain
        axes(grainIndex,1)=uniqueGrains(grainIndex);
        axes(grainIndex,2)= sqrt(length(XY)/3.1415) ;
        axes(grainIndex,3)=sqrt(length(XY)/3.1415);
        axes(grainIndex,4)=0;
        axes(grainIndex,5)=0;
        axes(grainIndex,6)=0;
        axes(grainIndex,7)=0;
        axes(grainIndex,8)=0;
        axes(grainIndex,9)=mu(1);
        axes(grainIndex,10)=mu(2);        
        continue;
    end

    princpAxes = [];
    princpAxes = R*[1 0;0 1];

    princpAxes(1,:)=princpAxes(1,:)*ellipse_t.a;
    princpAxes(2,:)=princpAxes(2,:)*ellipse_t.b;

    %hold on
    %quiver([mu(1);mu(1)],[mu(2);mu(2)],princpAxes(:,1),princpAxes(:,2),'k','LineWidth',1);
    
    axes(grainIndex,1)=uniqueGrains(grainIndex);
    axes(grainIndex,2)=ellipse_t.a;
    axes(grainIndex,3)=ellipse_t.b;
    axes(grainIndex,4)=ellipse_t.phi*360/(2*3.1415);
    axes(grainIndex,5)=princpAxes(1,1);
    axes(grainIndex,6)=princpAxes(2,1);
    axes(grainIndex,7)=princpAxes(1,2);
    axes(grainIndex,8)=princpAxes(2,2);
    axes(grainIndex,9)=mu(1);
    axes(grainIndex,10)=mu(2);
end
    
end

