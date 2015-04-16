function [dataLower] = lowerhemisphere(data)
% This function takes a set of unit vectors and converts them all to lower
% hemisphere. Data is input (and output) as an Nx3 matrix

x = data(:,1);
y = data(:,2);
z = data(:,3);

[tr,pl] = cart2sph(x,y,z);
[idx]   = find(pl(:) > 0);
pl(idx) = -pl(idx);
tr(idx) = tr(idx) + pi;

[x,y,z] = sph2cart(tr,pl,1);

dataLower = [x(:), y(:), z(:)];

end
