function [rotmatrix] = imrotate_3D(matrix,degrees)
%imrotate_3D Function to rotate a 3D image by a certain number of degrees

rotmatrix=zeros(size(matrix)); % midx and midy same for both

midx = round(size(rotmatrix,1)/2);
midy = midx;

for i=1:size(rotmatrix,1)
    for j=1:size(rotmatrix,2)

         x= (i-midx)*cosd(degrees)+(j-midy)*sind(degrees);
         y=-(i-midx)*sind(degrees)+(j-midy)*cosd(degrees);
         x=round(x)+midx;
         y=round(y)+midy;

         if (x>=1 && y>=1 && x<=size(matrix,2) && y<=size(matrix,1))
              rotmatrix(i,j,:)=matrix(x,y,:); % k degrees rotated image         
         end
    end
end

end