function [eulerAngles] = orientationmatrix2euler(g)
% Calculate Euler angles from a given rotation matrix using Bunge
% convention (passive, intrinsic ZXZ convention)
% euler angles are returned in radians

%% Check inputs
narginchk(1,1);

% g is given as Nx9 vector, see euler2rotationmatrix.m for details.
validateattributes(g,{'numeric'},{'ndims',2,'ncols',9});

%% Begin function

switch g(3,3)
    case -1 % Gimball lock, infinite solutions, set alpha = 0
        %... unclear how this affects crystal orientations
        %... Treat the singular case where beta = pi
        alpha = 0;
        beta = pi;
        gamma = atan2(R(1,2),R(1,1));
        warning('Gimball Lock encountered, alpha and gamma are not unique')
    case 1 % Gimball lock, infinite solutions, set alpha = 0,
        %... Treat the singular case where beta = 0
        alpha = 0;
        beta = 0;
        gamma = atan2(R(1,2),R(1,1));
        warning('Gimball Lock encountered, alpha and gamma are not unique')
    otherwise
        %... Treat the general case
        alpha = atan2(g(:,7),-g(:,8)); % atan2(g31,-g32)
        gamma = atan2(g(:,3), g(:,6)); % atan2(g13, g23)
        beta  = atan2(sqrt(g(:,3).^2 + g(:,6).^2),g(:,9)); % atan2(sqrt(g13^2 + g23^2), g33)
end

%... adjust angles from [-pi,pi] as output by atan2 to [0,2pi]
alpha = wrap(alpha,0,2*pi);
beta  = wrap(beta, 0,2*pi);
gamma = wrap(gamma,0,2*pi);


eulerAngles = [alpha beta gamma];

end

