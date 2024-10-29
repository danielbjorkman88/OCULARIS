function [out] = surfToObj(in)

[ii,~] = size(in.faces);

for jj = 1:ii
   out(jj).XYZ = ...
       [in.vertices(in.faces(jj,1),:);...
       in.vertices(in.faces(jj,2),:);...
       in.vertices(in.faces(jj,3),:)];
   
    
    out(jj).plane = fitPlane(out(jj).XYZ);
    
end



end