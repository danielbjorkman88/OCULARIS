function pointsOut = mapping(pointsIn,TMatrix)

% function for mapping of 3d points by roto-translation Matrix 
% pointsIn = Nx3 points in three dimensional coordinates
% TMatrix = 4x4 rt matrix

[n,p] = size(pointsIn);

pointsIn(:,4) = 1;

for i=1:n
   out(:,i) = TMatrix*pointsIn(i,:)';
end

pointsOut = out(1:3,:)';

end