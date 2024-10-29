function [angles] = getFixPntPA(point)
angles = [];


d2 = norm([point(1) point(2)]);
r2 = norm(point);

angles(1) = asind(d2/r2);
angles(2) = 2*atand(point(2)/(d2 + point(1)));

if angles(2)<0
    angles(2) = 360+angles(2);
end

if isnan(angles(2))
    angles(2) = 0;
end

end