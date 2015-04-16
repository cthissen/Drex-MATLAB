function [CrystalDirections,nDirections] = rotationmatrix2quartzdirections(R)
% this function converts euler angles into crystal directions appropriate
% for quartz

        % Plots c axes, <a> axes, and poles to r and m and ?
        nDirections = 5;
        
        % parse orientation matrix
        g11 = R(:,1); % R(:,1) being g21 seems to work for Lloyd samples...
        g12 = R(:,2);
        g13 = R(:,3);
        
        g21 = R(:,4);
        g22 = R(:,5);
        g23 = R(:,6);
        
        g31 = R(:,7);
        g32 = R(:,8);
        g33 = R(:,9);
        
        % get c axes
%         c = [g31, g32, g33];

        % get a axes by symmetry 
        a = sqrt(3)/2;        
        gp11 = -0.5*g11 + a*g21;
        gp21 = -0.5*g12 + a*g22;
        gp31 = -0.5*g13 + a*g23;
                
        gpp11= -0.5*g11 - a*g21;
        gpp21= -0.5*g12 - a*g22;
        gpp31= -0.5*g13 - a*g23;
        
        
        % get poles to m planes
        theta = [30, 90, 150, 210, 270, 330]*pi/180; % all m
        m = zeros([numel(g11),3,6]);
        for i = 1:numel(theta)
            cosTh = cos(theta(i));
            sinTh = sin(theta(i));
            m(:,1,i) = cosTh*g11 - sinTh*g21;
            m(:,2,i) = cosTh*g12 - sinTh*g22;
            m(:,3,i) = cosTh*g13 - sinTh*g23;
        end
            
        % get poles to r plane 
        for i = 1:6
%            r(:,1,i) = 0.5*(m(:,1,i) + g31);
%            r(:,2,i) = 0.5*(m(:,2,i) + g32);
%            r(:,3,i) = 0.5*(m(:,3,i) + g33);
           
           x = 0.5*(m(:,1,i) + g31);
           y = 0.5*(m(:,2,i) + g32);
           z = 0.5*(m(:,3,i) + g33);
           
           %... normalize to unit vector
           length = sqrt(sum(x.^2 + y.^2 + z.^2,2));
           r(:,1,i) = x./length;
           r(:,2,i) = y./length;
           r(:,3,i) = z./length;
            
        end

        % assign crystal axes
        CrystalDirections{1}.name        = 'c-axes';
        CrystalDirections{1}.nData       = numel(g11);
        CrystalDirections{1}.multiplicity= 1;        
        CrystalDirections{1}.unitVectors = [g31,g32,g33];  
        
        CrystalDirections{2}.name             = '<a> axes';
        CrystalDirections{2}.nData            = numel(g11);
        CrystalDirections{2}.multiplicity     = 3;                
        CrystalDirections{2}.unitVectors(:,1) = [g11; gp11; gpp11];
        CrystalDirections{2}.unitVectors(:,2) = [g21; gp21; gpp21];
        CrystalDirections{2}.unitVectors(:,3) = [g31; gp31; gpp31];

        CrystalDirections{3}.name             = 'poles to {m}'; % positive m only (1:3)
        CrystalDirections{3}.nData            = numel(g11);
        CrystalDirections{3}.multiplicity     = 3;
        CrystalDirections{3}.unitVectors(:,1) = [m(:,1,4); m(:,1,5); m(:,1,6)]; 
        CrystalDirections{3}.unitVectors(:,2) = [m(:,2,4); m(:,2,5); m(:,2,6)]; 
        CrystalDirections{3}.unitVectors(:,3) = [m(:,3,4); m(:,3,5); m(:,3,6)]; 
        
        CrystalDirections{4}.name             = 'poles to {z}';
        CrystalDirections{4}.nData            = numel(g11);
        CrystalDirections{4}.multiplicity     = 3;        
        CrystalDirections{4}.unitVectors(:,1) = [r(:,1,1); r(:,1,3); r(:,1,5)]; 
        CrystalDirections{4}.unitVectors(:,2) = [r(:,2,1); r(:,2,3); r(:,2,5)]; 
        CrystalDirections{4}.unitVectors(:,3) = [r(:,3,1); r(:,3,3); r(:,3,5)]; 
                
        CrystalDirections{5}.name             = 'poles to {r}';
        CrystalDirections{5}.nData            = numel(g11);
        CrystalDirections{5}.multiplicity     = 3;        
        CrystalDirections{5}.unitVectors(:,1) = [r(:,1,2); r(:,1,4); r(:,1,6)]; 
        CrystalDirections{5}.unitVectors(:,2) = [r(:,2,2); r(:,2,4); r(:,2,6)]; 
        CrystalDirections{5}.unitVectors(:,3) = [r(:,3,2); r(:,3,4); r(:,3,6)]; 
        

%... Bunge convention (ZXZ), passive rotation
% R(:,1)= c1.*c3 - c2.*s1.*s3;    %R11
% R(:,2)= c3.*s1 + c1.*c2.*s3;    %R12
% R(:,3)= s2.*s3;                 %R13
% 
% R(:,4)=-c1.*s3 - c2.*c3.*s1;    %R21
% R(:,5)= c1.*c2.*c3 - s1.*s3;    %R22
% R(:,6)= c3.*s2;                 %R23
% 
% R(:,7)= s1.*s2;                 %R31
% R(:,8)=-c1.*s2;                 %R32
% R(:,9)= c2;                     %R33


end