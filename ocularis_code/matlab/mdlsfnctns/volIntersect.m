function [out] = volIntersect(vol1,vol2,plot)
% vol1 = volMR vol2 = VolEP

out.faces = [];
out.vertices = [];
out.V = [];

vs1 = vol1.surf;
vs2 = vol2.surf;

v1 = surfToObj(vs1);
v2 = surfToObj(vs2);

[points12,points21] = inObject3d(v2,v1,2);

pt_inter = [points12;points21];

if isempty(pt_inter)
    out.V = 0;
    out.OL_MR = 0;
    return;
end

[K,V] = convhull(pt_inter(:,1),pt_inter(:,2),pt_inter(:,3));

out.faces = K;
out.vertices = pt_inter;


if V > vol1.V
    [K1,V1] = convhull(vs1.vertices(:,1),vs1.vertices(:,2),vs1.vertices(:,3));
    
    out.OL_MR = V/V1;
    out.V = V*(vol1.V/V1);
elseif V > vol2.V
    [K1,V2] = convhull(vs2.vertices(:,1),vs2.vertices(:,2),vs2.vertices(:,3));
    
    out.OL_MR = V/V2;
    out.V = V*(vol2.V/V2);
    
else    
    out.OL_MR = V/vol1.V;
    out.V = V;
end

% close
if plot
    figure, hold on, axis equal, grid on
    trisurf(vs1.faces,vs1.vertices(:,1),vs1.vertices(:,2),vs1.vertices(:,3),'FaceColor',[1 0 0]);
    trisurf(vs2.faces,vs2.vertices(:,1),vs2.vertices(:,2),vs2.vertices(:,3),'FaceColor',[0 0 1]);
    alpha(0.5)
%     plot3(points12(:,1),points12(:,2),points12(:,3),'+b')
%     plot3(points21(:,1),points21(:,2),points21(:,3),'+r')

%     trisurf(K,pt_inter(:,1),pt_inter(:,2),pt_inter(:,3))
    %     trisurf(K1,vol1.vertices(:,1),vol1.vertices(:,2),vol1.vertices(:,3))
end

out.V;

end