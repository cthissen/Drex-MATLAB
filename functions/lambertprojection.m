function [LambertProj] = lambertprojection(spherePoints,varargin)
% This function generates a set of approximately evenly spaced points on the 
% the upper (lower) hemisphere of a sphere.
% The function starts with equally spaced points in an equal area
% projection, and inverts the projection to get the 3D coordinates.
% 
% The input spherePoints^2 is the total number of points in the equal area
% projection. The actual points in the lower hemisphere of the stereonet
% will be somewhat lower. 
% 
% The output variables are formatted to spherePoints x spherePoints
% matrices
% 
% 
%... Lower Hemisphere (default)
%            %%%    %%%
%       %%%              %%%
% 
%   %%%                      %%%
% 
%  %%%                         %%%
% %               +z            +x
%  %%%                         %%%
% 
%  %%%                        %%%
% 
%     %%%                  %%%
% 
%           %%%  +y  %%%
% 
% %... Upper Hemisphere
%            %%%  +y  %%%
%       %%%              %%%
% 
%   %%%                      %%%
% 
%  %%%                         %%%
% %               +z            +x
%  %%%                         %%%
% 
%  %%%                        %%%
% 
%     %%%                  %%%
% 
%           %%%     %%%
%           

%% Check Inputs
narginchk(0,2);

validateattributes(spherePoints,{'numeric'},{'scalar','>',0})

%% Begin Function

% check for input
if isempty(varargin)
    % set default lower hemisphere
    hemisphere = 'lower';
else
    hemisphere = varargin{1};
end
    

%... build equal area projection
R = sqrt(2);
X = linspace(-R,R,spherePoints);
Y = X;
[X,Y] = meshgrid(X,Y); % X increases to right, 
Y = flipud(Y);         % Y increases up

%... calculate equivalent points in 3D
x = sqrt(1-(X.^2 + Y.^2)/4).*X;
switch hemisphere
    case 'lower'
        y = sqrt(1-(X.^2 + Y.^2)/4).* -Y;
        z = -(1 - (X.^2 + Y.^2)/2);

    case 'upper'
        y = sqrt(1-(X.^2 + Y.^2)/4).* Y;
        z = (1 - (X.^2 + Y.^2)/2);

    otherwise
        error('specify lower or upper for hemisphere');
end

% correct roundoff errors
x(1-(X.^2 + Y.^2)/4 < eps) = 0;
y(1-(X.^2 + Y.^2)/4 < eps) = 0;

% ensure unit vectors
mag = sqrt(x.^2 + y.^2 + z.^2);
x = x./mag;
y = y./mag;
z = z./mag;

% parse output
LambertProj.x = x;
LambertProj.y = y;
LambertProj.z = z;
LambertProj.X = X;
LambertProj.Y = Y;
LambertProj.R = R;

end

