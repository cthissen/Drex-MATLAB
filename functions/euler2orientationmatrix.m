function [g] = euler2orientationmatrix(eulerAngles)
% Euler2Rotation inputs a set of 3 Euler angles (alpha beta and gamma) and outputs
% a rotation matrix using the Bunge convention, (e.g. ZXZ, passive 
% intrinsic rotation).
%
% eulerAngles are input as radians

%% Check inputs
narginchk(1,1);

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

% parse inputs
Alpha = eulerAngles(:,1); % col vecs
Beta  = eulerAngles(:,2);
Gamma = eulerAngles(:,3);

%% Calculation
c1 = cos(Alpha);
s1 = sin(Alpha);

c2 = cos(Beta);
s2 = sin(Beta);

c3 = cos(Gamma);
s3 = sin(Gamma);

% p 34 Engler and Randle, 2010
%... Bunge convention (ZXZ), passive rotation
g(:,1)= c1.*c3 - c2.*s1.*s3;    %g11
g(:,2)= c3.*s1 + c1.*c2.*s3;    %g12
g(:,3)= s2.*s3;                 %g13

g(:,4)=-c1.*s3 - c2.*c3.*s1;    %g21
g(:,5)= c1.*c2.*c3 - s1.*s3;    %g22
g(:,6)= c3.*s2;                 %g23

g(:,7)= s1.*s2;                 %g31
g(:,8)=-c1.*s2;                 %g32
g(:,9)= c2;                     %g33


% % ... derivation
% syms c1 s1 c2 s2 c3 s3 real
% Ra = [c1 s1 0;
%      -s1 c1 0;
%        0  0 1];
%    
% Rb = [1 0 0 ;
%       0 c2 s2;
%      0 -s2 c2];
%  
% Rg = [c3 s3 0;
%      -s3 c3 0;
%        0  0 1];
%    
% R = Rg*Rb*Ra   

% Drexv2b.f90
%       acs0(i,1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
%       acs0(i,1,2)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
%       acs0(i,1,3)=SIN(phi2)*SIN(theta)
% 
%       acs0(i,2,1)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
%       acs0(i,2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
%       acs0(i,2,3)=COS(phi2)*SIN(theta)
% 
%       acs0(i,3,1)=SIN(theta)*SIN(phi1)
%       acs0(i,3,2)=-SIN(theta)*COS(phi1)
%       acs0(i,3,3)=COS(theta)

end
