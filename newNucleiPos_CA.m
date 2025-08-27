function [centers,centerIndex] = newNucleiPos_CA(grid,xgrid,ygrid,currentTemps,nucleationData,dt,dx)
    %Function to return the positions of the new nuclei given the current state fo the simulation    

    z=nucleationData.z;
    K1=nucleationData.K1;
    G_nucleation=nucleationData.G_nucleation;
    k=nucleationData.k;
    Tm=nucleationData.Tm;
    nucleationRateSTD=nucleationData.nucRateSTD;
    shapeFactor=nucleationData.shapeFactor;

    %Calc total number of new nuclei from homogenous nucleation
    %newNuclei=round(normrnd(nucleationRate,nucleationRateSTD)*dt);
    
    centers = [];
    centerIndex=[];

    Gstar = K1./(Tm-currentTemps).^2; %driving force for nucleation due to undercooling
    Nhomo=z.*exp(-Gstar./(currentTemps*k)).*exp(-G_nucleation./(currentTemps.*k)); %homogeneous nucleation rate

    gridZero= grid==0; %find all the spots open to nucleation
    probNucleation = Nhomo.*dt.*(dx^2).*gridZero; %probability of nucleation considering which spots of the grid are liquid
    
    %SETTING HOMOGENOUS NUCLEI TO ZERO
    %probNucleation = zeros(size(grid));

    %Perform heterogeneous nucleation
    hetPos=[]; %brute force check for heterogenous nucleation
    for row = 1:length(grid)
        for col = 1:length(grid)
            if grid(row,col)==0 %ensure it's unoccupied
                if nucleation_Boundary(col,row,grid)==1
                    probNucleation(row,col)=z*exp(-Gstar(row,col)*shapeFactor/(currentTemps(row,col)*k))*exp(-G_nucleation/(currentTemps(row,col)*k))*dt*(dx^2);
                    hetPos = [hetPos;row,col];
                end
            end
        end
    end
    randNum=rand(length(grid),length(grid)); %random numbers for each x,y location
    nucleiCheck=probNucleation>randNum; %find spots with sufficient probability of nuclei

    %Check if heterogenous nucleation occured
%     fullSize = size(hetPos);
%     if fullSize(1)>1
%         for i = 1:length(hetPos)
%             if nucleiCheck(hetPos(i,1),hetPos(i,2))==1
%                 fprintf("Het. nucleation at %i,%i\n",hetPos(i,1),hetPos(i,2));
%             end
%         end
%     elseif fullSize(1)==1
%         fprintf("Het. nucleation at %i,%i\n",hetPos(1),hetPos(2));
%     end
        
    openIndices = find(nucleiCheck==1);
    centerIndex = openIndices;
    centers = [xgrid(centerIndex),ygrid(centerIndex)];

end

