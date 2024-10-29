function drawAxis3datPosandOr(origin,V,varargin)
%DRAWAXIS3D Draw a coordinate system and an origin
%
%   drawAxis3d
%	Adds 3 cylinders to the current axis, corresponding to the directions
%	of the 3 basis vectors Ox, Oy and Oz.
%	Ox vector is red, Oy vector is green, and Oz vector is blue.
%
%   drawAxis3d(L, R)
%   Specifies the length L and the radius of the cylinders representing the
%   different axes.
%
%   Example
%   drawAxis
%
%   figure;
%   drawAxis(20, 1);
%
%   See also
%   drawAxisCube
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2007-08-14,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.

% geometrical data
origin = origin;
v1 = V(1,:);
v2 = V(2,:);
v3 = V(3,:);

% default parameters
L = 1;
r = L/10;

% extract parameters from input
if ~isempty(varargin)
	L = varargin{1};
end
if length(varargin)>1
	r = varargin{2};
end

% draw 3 cylinders and a ball
hold on;
drawCylinder([origin+v1*1 origin+v1*L r], 16, 'facecolor', [.8 .8 .8], 'edgecolor', 'black');
drawCylinder([origin+v2*1 origin+v2*L r], 16, 'facecolor', [.8 .8 .8], 'edgecolor', 'black');
drawCylinder([origin+v3*1 origin+v3*L r], 16, 'facecolor', [.8 .8 .8], 'edgecolor', 'black');
% drawSphere([origin 2*r], 'faceColor', 'black');

