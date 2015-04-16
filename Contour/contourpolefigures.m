function [hFig] = contourpolefigures(eulerAngles,mineral,method,varargin)
% For an input of euler angles (Bunge convention) and a mineral name
% ('olivine' or 'quartz'), this function plots the pole figures using Kamb
% or Gaussian contouring. The function is currently setup only for
% non-polar data, and calculates the density in a single hemisphere.
%
% Input parameters: eulerAngles is an Nx3 matrix of Bunge convention Euler
% angles in radians.
%
% mineral is 'olivine' or 'quartz' and determines which crystallographic
% axes are plotted. Olivine pole figures include [100], [010], and [001]
% directions. Quartz pole figures include c-axes, a-axes, and poles to m,r,
% and z planes.
%
% method is the contouring method, and is either Kamb or Guassian
% 
% The figure handle can be input as varargin. Otherwise a new figure is
% created. Addition variable input arguments are the stereogram hemisphere,
% the number of points to use in each dimension of the stereogram, and
% options for the colormap.
%
%
% Version 1.0. CJT Jan 27, 2015. First working version. 
% Version 1.01 CJT Feb  4, 2015. Included optional figure handle input
% Version 1.1  CJT Mar 18, 2015. Added gaussian smoothing
%
%... default coordinate frame
%            %%%    %%%
%       %%%              %%%
% 
%   %%%                      %%%
% 
%  %%%                         %%%
% %               +z             +x
%  %%%                         %%%
% 
%  %%%                        %%%
% 
%     %%%                  %%%
% 
%           %%%  +y  %%%

% Dependencies:
% addpath('/Volumes/Research/MATLAB/');
% addpath('/Volumes/Research/MATLAB/EulerAngleFunctions/');
% addpath('/Volumes/Research/MATLAB/Stereogram/');
% addpath('/Volumes/Research/MATLAB/cmapscale/');
% addpath('/Volumes/Research/MATLAB/colormaps/');
% addpath('/Volumes/Research/MATLAB/dispstat/');
% dispstat('','init'); % One time only initialization 


%% Check inputs
narginchk(2,7); % requires 2 argument, with 3 additional optional arguments

% check eulerAngles is Nx3
[~,nDim] = size(eulerAngles);
if nDim ~= 3
    error('eulerAngles must be input as Nx3 matrix');
end

% check range of eulerAngles
maxE = max(eulerAngles,[],1);
minE = min(eulerAngles,[],1);
if any(maxE > 2*pi) || any(minE < 0) || maxE(2) > pi
    error('Expected euler angles be in range 0 < eulerAngles(:,1), eulerAngles(:,3) < 2pi and 0 < eulerAngles(:,2) < pi');
end

% validate method of contouring input
validatestring(method,{'Kamb','Gaussian'});

% Parse input and define default values
%... colormap options for cmapscale.m
ColorOpts.colorRamp0 = hklcolor(256);
ColorOpts.nTics = 3;
ColorOpts.centerVal = [];
ColorOpts.factor = 0.8;

optArgs = {1,'lower',101,ColorOpts};
nArgsIn = find(~cellfun(@isempty,varargin));


%...copy cells from varargin to optArgs
optArgs(nArgsIn) = varargin(nArgsIn);

[f,hemisphere,nSpherePoints,ColorOpts] = optArgs{:};


%% Start Calculation
% Convert euler angles to unit vectors
g = euler2orientationmatrix(eulerAngles); % g is arranged as nx9 vector. 

% convert dir cosines to relevant crystal directions
[CrystalDirections,nDirections] = rotationmatrix2crystaldirections(g,mineral);

% Calculate points on the sphere
[SphereProj] = lambertprojection(nSpherePoints,hemisphere);



% operate on each crystal direction
for i = 1:nDirections
    % convert axes to lower hemisphere (find antipodes)
%     CrystalDirections{i}.unitVectors = lowerhemisphere(CrystalDirections{i}.unitVectors);

    % Get Kamb contours levels
    switch method
        case 'Kamb'
            CrystalDirections(i) = kambcontour(CrystalDirections(i),SphereProj);
        
        
        case 'Gaussian'
            CrystalDirections(i) = gaussiancontour(CrystalDirections(i),SphereProj);

    end
 
end 

% Make figure
hFig = figure(f); clf
hFig = makepolefigures(CrystalDirections,mineral,SphereProj,hFig,ColorOpts);

%% LINT: This will plot all datapoints in 3D
% figure(3); clf
% subplot(1,3,1)
% scatter3(CrystalDirections{1}.unitVectors(:,1),CrystalDirections{1}.unitVectors(:,2),CrystalDirections{1}.unitVectors(:,3))
% xlabel('x')
% axis equal
% 
% subplot(1,3,2)
% scatter3(CrystalDirections{2}.unitVectors(:,1),CrystalDirections{2}.unitVectors(:,2),CrystalDirections{2}.unitVectors(:,3))
% xlabel('x')
% axis equal
% 
% subplot(1,3,3)
% scatter3(CrystalDirections{3}.unitVectors(:,1),CrystalDirections{3}.unitVectors(:,2),CrystalDirections{3}.unitVectors(:,3))
% xlabel('x')
% axis equal

% scatter3(CrystalDirections.data{1}.x(:),CrystalDirections.data{1}.y(:),CrystalDirections.data{1}.z(:))

%% 
% keyboard

end




        