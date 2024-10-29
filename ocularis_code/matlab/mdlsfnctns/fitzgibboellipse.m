function A = fitzgibboellipse(XY)
A = [];
centroid =[0,0];
D  = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
      (XY(:,2)-centroid(2)).^2,XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];

S = D'*D;
C(6,6) = 0;C(1,3) = 2;C(2,2) = -1;C(3,1) = 2;
cc = (S)\C;
if or(any(any(isnan(cc))),any(any(isinf(cc))))
 
    clearvars -except A
    return
end
[gevec,geval] = eig(cc);
[posR,posC] = find(geval>0 & ~isinf(geval));
    
if isempty(posC)
    clearvars -except A;
    return;
else
    A = gevec(:,posC(1));
    clearvars -except A
end
end