function [diode3D] = getFixPnt3D(pol,azim)

% get 3D position of fixation diode

d__ = 132.5;
rPolar = d__*tand(pol);
diode3D = [rPolar*cosd(azim) rPolar*sind(azim) d__];


end