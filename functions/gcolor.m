function h = gcolor(varargin)
%GCOLOR Pseudocolor (checkerboard) plot for grid-registered data.
%   GCOLOR(C) is a pseudocolor or "checkerboard" plot of matrix C,
%   where the values in C are centered on the grid nodes.
%   The values of the elements of C specify the color in each
%   cell of the plot. In the default shading mode, 'faceted',
%   each cell has a constant color and the last row and column of
%   C are not used. With shading('interp'), each cell has color
%   resulting from bilinear interpolation of the color at its 
%   four vertices and all elements of C are used. 
%   The smallest and largest elements of C are assigned the first and
%   last colors given in the color table; colors for the remainder of the 
%   elements in C are determined by table-lookup within the remainder of 
%   the color table.
%
%   GCOLOR(X,Y,C), where X and Y are vectors or matrices, makes a
%   pseudocolor plot on the grid defined by X and Y. Values in 
%   C correspond to the grid nodes indicated by X and Y.
%
%   GCOLOR(AX,..) plots into AX instead of GCA.
%
%   H = GCOLOR(...) returns a handle to a SURFACE object.
%
%   GCOLOR is really a SURF with its view set to directly above.
%
%   See also PCOLOR CAXIS, SURF, MESH, IMAGE, SHADING.

%-------------------------------
%   Additional details:
%
%   GCOLOR sets the View property of the SURFACE object to directly 
%   overhead.
%
%   If the NextPlot axis property is REPLACE (HOLD is off), GCOLOR resets 
%   all axis properties, except Position, to their default values
%   and deletes all axis children (line, patch, surf, image, and 
%   text objects).  View is set to [0 90].

%   Adapted from the Matlab GCOLOR program, which is covered
%   by the following copyright and revision number.
%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 5.9.4.6 $  $Date: 2011/07/25 03:49:30 $
%   J.N. Little 1-5-92

%% Start function
% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,3,nargs,'struct'))

% do error checking before calling newplot. This argument checking should
% match the surface(x,y,z) or surface(z) argument checking.
if nargs == 2
  error(message('gcolor:InvalidNumberOfInputs'))
end
if isvector(args{end})
  error(message('gcolor:NonMatrixColorInput'));
end
if nargs == 3 && LdimMismatch(args{1:3})
  error(message('gcolor:InputSizeMismatch'));
end
for k = 1:nargs
  if ~isreal(args{k})
    error(message('gcolor:NonRealInputs'));
  end
end

cax = newplot(cax);
hold_state = ishold(cax);

if nargs == 1
    x = args{1};
    [m,n] = size(x);
    hh = surface(zeros(m+1,n+1),...
        padarray(c,[1 1],'replicate','post'), ...
        'parent',cax);
    lims = [ 1 n 1 m];
elseif nargs == 3
    [x,y,c] = deal(args{1:3});
    lims = [min(x(:)),max(x(:)),min(y(:)),max(y(:))];
    [m,n] = size(c);
    dX = (x(1,end) - x(1,1))/(n-1);
    dY = (y(end,1) - y(1,1))/(m-1);
    xNew = x(1,1) - dX/2 : dX : x(1,end) + dX/2;
    yNew = (y(1,1) - dY/2 : dY : y(end,1) + dY/2)';
    hh = surface(xNew,yNew,zeros(m+1,n+1), ...
        padarray(c,[1 1],'replicate','post'), ...
        'parent',cax);
end
if ~hold_state
    set(cax,'View',[0 90]);
    set(cax,'Box','on');
    axis(cax,lims);
end
if nargout == 1
    h = hh;
end

function ok = LdimMismatch(x,y,z)
[xm,xn] = size(x);
[ym,yn] = size(y);
[zm,zn] = size(z);
ok = (xm == 1 && xn ~= zn) || ...
     (xn == 1 && xm ~= zn) || ...
     (xm ~= 1 && xn ~= 1 && (xm ~= zm || xn ~= zn)) || ...
     (ym == 1 && yn ~= zm) || ...
     (yn == 1 && ym ~= zm) || ...
     (ym ~= 1 && yn ~= 1 && (ym ~= zm || yn ~= zn));
