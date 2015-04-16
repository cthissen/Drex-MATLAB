function [colorRamp1,tics] = cmapscale(z,colorRamp0,varargin)
% cmapscale does data-driven rescale of the color-ramp distribution
%    COLORRAMP1 = cmapscale(Z, COLORRAMP0) rescales COLORRAMP0 so that it matches
%    the cumulative distribution for the data Z.
%
%    COLORRAMP1 = cmapscale(Z,COLORRAMP0,FACTOR) adjusts the contrast from
%    linear mapping (no change) when FACTOR = 0, to uniform mapping (maximum
%    contrast) when FACTOR = 1. The default is uniform mapping, which is
%    equivalent to histogram equalization.
%
%    COLORRAMP1 = cmapscale(Z,COLORRAMP0,FACTOR,Z0) rescales so that the
%    value Z0 lies at the center of the color ramp (e.g. to center the color
%    ramp around zero). Default is no centering.
%
%    [COLORRAMP1,TICS] = cmapscale(Z,COLORRAMP0,FACTOR,Z0,NTICS) returns an array
%    TICS, with size nTics x 2. The colorbar is more informative if it uses the
%    original color ramp. TICS are used to locate tics and tic labels relative
%    to COLORRAMP0. The default is no tics. The first column of TICS
%    has coordinates between 0 and 1 for the initial color ramp, and
%    the second column contains the respective Z values.
% 
%    By Mark Brandon, Yale University, August, 2012.
%    Inspired by Eisemann, M., Albuquerque, G., and Magnor, M., 2011, 
%       Data driven color mapping: Proceedings EuroVA.
% 
%    1-29-15. CJT. Modified calculation of tics to allow extrapolation (and
%    remove NaN by adding 'linear','extrap' to cLTics = interp1().
% 
%    3-13-15. CJT. Added a routine to check for two or fewer unique values in the
%    data Z and modify colorRamp1 and tics accordingly.
%
%    3-25-15. CJT. Added an addition varargin that allows the user to
%    set the limits on the colormap. For example, [0,1] sets the colormap
%    limits to 0 and 1. The default is empty, which sets the color limits
%    based on the data. 
%% Define default values and validate input arguments
narginchk(2,6); % at least 2 but no more than 6 input arguments

% Define default values
factor = 1; % uniform mapping (maximum contrast)
z0    = []; % no centering value
nTics = []; % no set number of tics
cLims = []; % let color limits be defined by data
optargs = {factor, z0, nTics,cLims};

numvarargs = length(varargin); % number of variable input arguments
optargs(1:numvarargs) = varargin; % transfer varargin into opt args
[factor, z0, nTics, cLims] = optargs{:}; % parse varagin in into variables used later


% Validate input arguments
%... Check for contrast factor
if factor<0 || factor>1
    error('Factor must be in the range 0 to 1.');
end

%... Check for minimum number of unique z values
nUnique = numel(unique(z(~isnan(z))));
if nUnique <= 2
    fprintf('cmapscale.m warning: Data values have < 3 unique values.\n');
end

%... Check that centering value is within data range (or within cLims)
if isempty(cLims)
    % check if z0 is between data
    if ~isempty(z0) && (z0 <= min(z(:)) || z0 >= max(z(:)))
        z0 = [];
        warning('Centering value z0 is outside range for data values in z. Resetting z0=[]');
    end
else
    % check if z0 is between color limits
    if ~isempty(z0) && (z0 <= cLims(1) || z0 >= cLims(2))
        error('Centering value z0 is outside range for color limits.');
    end    
end

%... Check that nTics is a positive number
validateattributes(nTics,{'numeric'},{'positive'},'cmapscale','nTics')

%... Check that cLims make sense
if ~isempty(cLims)
    validateattributes(cLims,{'numeric'},{'increasing','size',[1,2]},'cmapscale','cLims');
end
%% Start Function
%... Recast contrast factor as a stretching parameter s, and limit
% its range to avoid numerical problems.
s = tan(pi*factor/2);
if s>1e4, 
    s = 1e4;
end

%... Sort and normalize z to a linear scale cL, and a uniform scale cU
if ~isempty(cLims)
    % set limits for colormap based on cmap inputs    
    z = z(:);
    z = sort(z(~isnan(z(:))));
    zMin = cLims(1);
    zMax = cLims(2);
    
else
    % set limits for colormap based on data
    if nUnique > 1    
        z = z(:);
        z = sort(z(~isnan(z(:))));
        zMin = z(1);
        zMax = z(end);
    else
        z = z(~isnan(z(:)));
        zMin = mean(z(~isnan(z(:))))-1e-7;
        zMax = zMin + 2e-7;
    end
end

%... Limit the size of the distribution to 1000 (to avoid memory problems)
if length(z) > 1000
    z = z(round((0:1/999:1)*(length(z)-1))+1);
end
cL = (z-zMin)/(zMax-zMin);
m=size(cL,1);
cU=(0:m-1)'./(m-1);

%... Rounding to 1e-7, which ensures dynamic range but avoids round off problems
cL = round(cL*1e7)*1e-7;

%... Create a new color scale cS using weighted sum of cL and cU
if isempty(z0)
    %... Color scale with no specified center point
    cS = (cL + s^2*cU)./(1 + s^2);
else
    %... Color scale with specified center point 
    i0 = sum(z(:)<=z0);
    if i0 == numel(z);
        warning('cmapscale.m: defined center value is greater than all data values');
        i0 = i0-1;
    end
    if z(i0)==z0
        %... Center point matches existing value in z
        cL0 = cL(i0);
        cU0 = cU(i0);
    else
        %... Center point lies between values in z
        cL0 = (z0-zMin)/(zMax-zMin);
        cL = [cL(1:i0); cL0; cL(i0+1:m)];
        cU0 = cU(i0) + (cU(i0+1) - cU(i0))/2;
        cU = [cU(1:i0); cU0; cU(i0+1:m)];
        i0 = i0 + 1;
        m = m + 1;
    end
    %... Project the lower and upper halves of the data distributions onto the
    % new scale c1. This is done in two steps to maintain z0 at the center.
    % Project lower half
    cS = zeros(m,1);
    cS(1:i0) = 0.5*(cL(1:i0)/cL0 + s^2*cU(1:i0)/cU0)/(1 + s^2);
    cS0 = cS(i0);
    % Project upper half
    cS(i0:m) = cS0 + (1-cS0)* ...
        (cL(i0:m)/(1-cL0) + s^2*(cU(i0:m)-cU0)/(1-cU0))/(1 + s^2);
end


%... Remove duplicates in the z sequence (and cL sequence as well).
% Each set of duplicates is reduced to a single value in the cL sequence,
% and the average value in the cS sequence. The result is a strictly
% monotonic sequence, which is required when using interp1 for 
% interpolation.
% Create the adjaceny matrix A for z
A = true(m,m);
for i = 1:m
    A(i,:) = (cL==cL(i));
end
% Calculate means for duplicate sequences
w = sum(A,2);
cL = (A*cL)./w;
cS = (A*cS)./w;
clear A
% Consolidate into a unique sequence
[cL,ia,~] = unique(cL);
cS = cS(ia);
m=size(cL,1);

if nUnique > 1

    % Normalize to account for round-off errors
    cL = cL/cL(m);
    cS = cS/cS(m);
    %... Interpolate new colors for color ramp
    mRamp = size(colorRamp0,1);
    cLRamp = (0:mRamp-1)'/(mRamp-1);
    cSRamp = interp1(cL,cS,cLRamp);
    colorRamp1 = interp1(cLRamp,colorRamp0,cSRamp);
    colorRamp1 = round(colorRamp1.*1e7)*1e-7;
    
    %... Calculate tick values
    if ~isempty(nTics) && nTics~=0
        cSTics = (0:nTics-1)'/(nTics-1);
        cLTics = interp1(cS,cL,cSTics,'linear','extrap');
        zTics = zMin + cLTics*(zMax-zMin);
        tics = [cSTics,zTics];
    end
else
    % we don't want to normalize if there is only 1 value
    mRamp = size(colorRamp0,1);
    cLRamp = (0:mRamp-1)'/(mRamp-1);    
    colorRamp1 = interp1(cLRamp,colorRamp0,cS);
    
    if ~isempty(nTics) && nTics~=0
        cSTics = (0:nTics-1)'/(nTics-1);
        cLTics = cL*ones(nTics,1);
         zTics = zMin + cLTics*(zMax-zMin);
        tics = [cSTics,zTics];

    end
end
    



end
