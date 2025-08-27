function [nuclei] = growNuclei(nuclei,growthRate,dt)
%growNuclei Function to grow all of the nuclei

% column 1-3: x,y,z position of the nuclei
% column 4: index of the cell that the nuclei belongs to
% column 5-7: euler angles of the nuclei
% column 8: initial nucleation time of the nuclei
% column 9: size of the nuclei

if size(nuclei,1)>0 %only perform if there are any nuclei to grow
    nuclei(:,9)= nuclei(:,9) + growthRate*dt;
end

end