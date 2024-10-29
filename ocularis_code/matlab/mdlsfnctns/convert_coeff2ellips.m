function [elli] = convert_coeff2ellips(coeff)

a = coeff(1);
b = coeff(2);
c = coeff(3);
d = coeff(4);
e = coeff(5);
f = coeff(6);


ac4_b2 = 4*a*c-b^2;


q =(64*f*(ac4_b2)- 64*(a*e^2-b*d*e+c*d^2))/...
    (ac4_b2^2);

s = 1/4*(sqrt(abs(q)*sqrt(b^2+(a-c)^2)));

r_max = 1/8*sqrt(2*abs(q)* sqrt(b^2+(a-c)^2) - 2*q*(a+c));

% r_max2 = sqrt(2*((a*)))


r_min = sqrt(r_max^2-s^2);

x_c = (b*e - 2*c*d)/(4*a*c - b^2);
y_c = (b*d - 2*a*e)/(4*a*c - b^2);

theta = 0;
if((q*a-q*c)==0 && q*b==0)
    theta = 0;
elseif((q*a-q*c)==0 && q*b>0)
    theta = rad2deg(pi/4);
elseif((q*a-q*c)==0 && q*b<0)
    theta = rad2deg(pi*3/4);
elseif((q*a-q*c)>0 && q*b>=0)
    theta = 0.5*atand(b/(a-c));
elseif((q*a-q*c)>0 && q*b<0)
    theta = 0.5*atand(b/(a-c))+180;
elseif((q*a-q*c)<0)
    theta = 0.5*atand(b/(a-c))+90;
end

elli = [x_c,y_c,r_max,r_min,theta];
clearvars -except elli
end