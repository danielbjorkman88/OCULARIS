function abs_volume = Vol_Polygon_pyramids(Polygon)
% calculate the volume of an triangled polygon
% 1) create center of gravity
% 2) create normal to the facets == distance t
% 3) basic area is created as cross product of facets
% 4) Vol = 1/3 basic area * height(t)


TR_Ell = Polygon;
% tic
E1 = TR_Ell.Points(TR_Ell.ConnectivityList(:,2),:) -TR_Ell.Points(TR_Ell.ConnectivityList(:,1),:);
E2 = TR_Ell.Points(TR_Ell.ConnectivityList(:,3),:) -TR_Ell.Points(TR_Ell.ConnectivityList(:,1),:);

N1 = cross(E1,E2,2);
N1_norm = N1./repmat(sqrt(dot(N1,N1,2)),1,3);
d1 = -dot(N1_norm,TR_Ell.Points(TR_Ell.ConnectivityList(:,1),:),2);

C = repmat(sum(TR_Ell.Points(:,:),1)./size(TR_Ell.Points,1),size(N1_norm,1),1);

t_mesh = -(d1+ dot(N1_norm,C,2))./dot(N1_norm,N1_norm,2);

basic_area = sqrt(dot(N1,N1,2))./2; % fuer die Fläche braut man N1! sonst norm
abs_volume = sum(abs(basic_area.*t_mesh./3));
% toc
% % abs() für t ????
% r1 = 11.6781; 
% r2 = 10.0543;
% 
% vgl_vol = 4/3*pi*r1*r1*r2


% tic
% %%%% das geht auch kürzer: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TR_Ell = Polygon;
% v0 = TR_Ell.Points(TR_Ell.ConnectivityList(:,1),:);
% v1 = TR_Ell.Points(TR_Ell.ConnectivityList(:,2),:);
% v2 = TR_Ell.Points(TR_Ell.ConnectivityList(:,3),:);
% sum(dot(v0, cross(v1, v2,2),2)./6)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toc