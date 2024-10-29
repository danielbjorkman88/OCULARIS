function [int_points12, int_points21] = inObject3d (Obj1,Obj2,f)
%% Determina i punti di una struttura interni all'altra e viceversa
% per f = 1 effettua il controllo solo dei punti Obj2 dentro Obj1
% per f = 2 effettua il controllo incrociato

int_points12 = [];
int_points21 = [];


aa = [];
for jj = 1:size(Obj1,2)
    aa = cat(1,aa,Obj1(jj).XYZ);
end
% aa = unique(aa,'rows');

for kk = 1:size(aa,1)
    n = isnan(aa(kk,:));
    if sum(n)~=0
        aa(kk:end,:) = [];
        break
    end
end

bb = [];
for jj = 1:size(Obj2,2)
    bb = cat(1,bb,Obj2(jj).XYZ);
end

% bb = unique(bb,'rows');

if (f == 2)
    [A1,b1] = vert2lcon (aa);
    A1 = ceil(A1*10000)/10000;
    b1 = ceil(b1*10000)/10000;
    
    int_1 = all(A1*(bb')<=b1);
    int_points12 (:,1) = bb(find(int_1),1);
    int_points12 (:,2) = bb(find(int_1),2);
    int_points12 (:,3) = bb(find(int_1),3);
    
    [A2,b2] = vert2lcon (bb);
    A2 = ceil(A2*10000)/10000;
    b2 = ceil(b2*10000)/10000;
    
    int_2 = all(A2*(aa')<=b2);
    int_points21 (:,1) = aa(find(int_2),1);
    int_points21 (:,2) = aa(find(int_2),2);
    int_points21 (:,3) = aa(find(int_2),3);
    
elseif f == 1
    
    [A1,b1] = vert2lcon (aa);
    A1 = ceil(A1*10000)/10000;
    b1 = ceil(b1*10000)/10000;
    
    int_1 = all(A1*(bb')<=b1);
    int_points12 (:,1) = bb(find(int_1),1);
    int_points12 (:,2) = bb(find(int_1),2);
    int_points12 (:,3) = bb(find(int_1),3);
    
end
end