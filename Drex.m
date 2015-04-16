function [] = DRex()
% This code simulates the development of LPO textures in a deforming
% olivine aggregate. This is a modified version of the popular fortran code
% D-Rex, originally released by Ed Kaminski, Neil Ribe, and others. 
% 
% The code has been modified from the original fortran version. This
% version of the code only tracks the development of crystal sizes and
% orientations. All path-line integration and seismic anistropy
% calculations present in the original fortran version have been removed.
% All calculations relating to enstatite have also been removed. The
% code is suitable for simulating olivine fabrics, and can easily be
% modified for calculating fabric development for a given velocity field.
% 
% Other modifications include:
% 1. The initially random LPO is now generated using a random stream, and
% follows the method in Morawiec, A. (2003). Orientations and rotations.
% Springer-Verlag.
% 2. The code warns the user when the rotation matrix used to update a
% crystal orientation is not a proper rotation matrix. 
% 3. An indexing error in the original fortran code related to the
% calculation of dislocation densities has be fixed. 
% 
% The theoretical development of the model is described in the following
% publications: 
% Kaminski, E., & Ribe, N. M. (2001). A kinematic model for 
% recrystallization and texture development in olivine polycrystals.
% Earth and Planetary Science Letters, 189(3), 253-267.
% 
% Kaminski, E., Ribe, N. M., & Browaeys, J. T. (2004). D-Rex, a program 
% for calculation of seismic anisotropy due to crystal lattice preferred 
% orientation in the convective upper mantle. Geophysical Journal
% International, 158(2), 744-752.
%
% %%%%%%%%%%%%%
%     USAGE
% %%%%%%%%%%%%%
% INPUT: 
% The code does not require input arguments. The user specifies two types
% of parameters in the following sections. The first type of parameters
% relates to the olivine aggregate, such as the number of olivine crystals
% and the stress exponent. These parameters are defined in the next section
% as part of the "Grain" structure. Default values and alternative options
% are also discussed.
% 
% The second type of input parameter relates to the imposed deformation.
% The user may choose from among the pre-defined deformation gradient
% tensors by setting the "deformationSymmetry" variable to one of the
% following strings: axisymmetricCompression, axisymmetricExtension,
% orthorhombic, simpleShear, triclinic, or noDeformation. Custom
% deformation gradient tensors can also be added in the section labelled
% "Define Deformation Gradient Tensor". The number of integration steps
% must also be specified.
% 
% OUTPUT:
% The code outputs several text files and figures. Five text files are
% created. The first is an info file that lists important parameters and
% final finite strain results for the aggregtate. The second text file
% lists the final orientations (ZXZ, Bunge convention, in radians),
% and the volume fraction of each crystal. Three other text files are
% also output that contain only euler angles. The unweighted file includes
% a single measurement for all crystals that have volume fraction greater
% than zero. This is akin to making a single measurement of each grain. The
% volume weighted file includes repeated measurements for the same grain scaled by the final volume fraction. 
% This is akin to making measurements on a predefined grid. The inverse
% volume weighted file includes repreate measurements for each grain
% 
% The code also outputs a figure containing pole figures for the
% [100],[010], and [001] axes. The default coordinate frame for the pole
% figures is defined as follows:
%            %%%     %%%
%       %%%              %%%
% 
%   %%%                      %%%
% 
%  %%%                         %%%
% %              +Z             +X % 
%  %%%                         %%%
% 
%  %%%                        %%%
% 
%     %%%                  %%%
% 
%           %%%  +Y  %%%
% The pole figures are saved using the export_fig toolbox. A version of this
% toolbox is included.

% V 0.4c. 4-7-15 CJT. First Working Version.
close all; clear all; clc

% Dependencies: contourpolefigures.m V1.0
addpath('Contour/');
addpath('export_fig/');
addpath('dispstat/');
addpath('colormaps/');
addpath('functions/');
version = '0.4c';

% keyboard
%% User Controlled Input Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
% Grain input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
Grain.nGrains   = 500;              % number of olivine crystals
Grain.pctOli    = 1;                % fraction of aggregate that is olivine
Grain.tau       = [1,2,3,1e60];     % CRSS (relative) for slip systems
Grain.mob       = 125;              % grain mobility parameter (125 is recommended by Kaminski et al 2004)
Grain.chi       = 0;                % 0.3 threshold volume fraction for activation of grain boundary sliding (Kaminski et al, 2004)
Grain.lambda    = 5;                % 5 nucleation parameter
Grain.stressExp = 3.5;              % stress exponent

%%%%%%%%%%%%%%%%%%%%%%%%
% Flow input parameters
%%%%%%%%%%%%%%%%%%%%%%%%
Flow.tSteps = 10;
Flow.deformationSymmetry = 'axisymmetricCompression';

% Various olivine "types" may simulated by changing the Grain.tau variable
% to the appropriate relative critical resolved shear stress:
% A-type: [1,2,3,1e60];
% B-type: [3,2,1,1e60];
% C-type: [3,2,1e60,1];
% D-type: [1,1,3,1e60]; 
% E-type: [3,1,2,1e60];
% Additional information on olivine types can be found in: 
% Karato, S. I., Jung, H., Katayama, I., & Skemer, P. (2008). Geodynamic
% significance of seismic anisotropy of the upper mantle: new insights from
% laboratory studies. Annu. Rev. Earth Planet. Sci., 36, 59-95.

% Flow.deformationSymmetry must be set to one of the following strings:
% axisymmetricCompression, axisymmetricExtension, orthorhombic, 
% simpleShear, triclinic, or noDeformation
%% Check Grain and Flow parameters
validateattributes(Grain,          {'struct'},{});
validateattributes(Grain.nGrains,  {'numeric'},{'scalar','integer','nonnegative'});
validateattributes(Grain.pctOli,   {'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(Grain.tau,      {'numeric'},{'nrows',1,'ncols',4});
validateattributes(Grain.mob,      {'numeric'},{'scalar','nonnegative'});
validateattributes(Grain.chi,      {'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(Grain.lambda,   {'numeric'},{'scalar','nonnegative'});
validateattributes(Grain.stressExp,{'numeric'},{'scalar','nonnegative'});

validateattributes(Flow,                     {'struct'},{});
validateattributes(Flow.tSteps,              {'numeric'},{'scalar','integer','nonnegative'});
validateattributes(Flow.deformationSymmetry, {'char'},{});

%% Define Deformation Gradient Tensor

switch Flow.deformationSymmetry
    
    case 'axisymmetricCompression' % Wk=0
        Sz = 0.42^2;
        Sx = sqrt(1/Sz);
        Sy = Sx;
        
        F = [sqrt(Sx),        0,       0;
                    0, sqrt(Sz),       0;
                    0         0, sqrt(Sy)];
              
    case 'axisymmetricExtension' % Wk=0        
        Sz = 1/0.42;
        Sx = sqrt(1/Sz);
        Sy = Sx;
        
        F = [sqrt(Sx),        0,       0;
                    0, sqrt(Sz),       0;
                    0         0, sqrt(Sy)];

    case 'orthorhombic' % Wk=0        
        Sx = 2;
        Sy = 1;
        Sz = 0.5;

        F = [sqrt(Sx),        0,       0;
                    0, sqrt(Sy),       0;
                    0         0, sqrt(Sz)];        
               
    case 'simpleShear' % Wk =/= 0
        % shear direction is x, shear plane is xy
        
        gammaDeg = 50;
        F = [1, tand(gammaDeg), 0;
             0,              1, 0;
             0               0, 1];          
%         F = [1, 1.2, 0;
%              0,              1, 0;
%              0               0, 1];          
               
    case 'triclinicShear' % Wk =/= 0
        % shear direction is x, shear plane is xy
        gammaHzDeg = 50;
        gammaVtDeg = 50;
        shortening = sqrt(0.75);        
        
        % shear in the x-direction normal to y
        FHzShear = [1, tand(gammaHzDeg), 0;
                    0,              1, 0;
                    0               0, 1];  
        
        % shear in the z-direction, normal to y        
        FVtShear = [1, 0, 0;
                    0, 1, 0;
                    0  tand(gammaVtDeg), 1];  
                
        FCoax = [1, 0, 0;        
                 0, shortening, 0;
                 0, 0, 1/shortening]; 

        F = FHzShear * FVtShear * FCoax;
        
        
        % This type of deformation gradient tensor gives shear vertically,
        % horiztonally, with compression. Constant volume is maintained by
        % allowing the material to slip along the boundaries (see Davis and Titus,
        % 2011, Figure 6)
        %         v1 = F(1,2)
        %         v2 = F(2,2)-1
        %         v3 = 2*F(3,2)*((1+v2)/(2+v2))
        %         F = [1,   v1, 0;
        %              0,  1+v2, 0;
        %              0, (v3/2)*(2+v2)/(1+v2), 1/(1+v2)];
      
    case 'noDeformation'
        % generate uniform random distribution of euler angles
        [eulerAngles] = generaterandomLPO(Grain);
            
        fileName = ['drex_',deformationSymmetry,'.txt'];
        fid = fopen(fileName,'w');
        for i = 1:Grain.nGrains
            fprintf(fid,'%f %f %f\n', eulerAngles(i,:)*180/pi);
        end
        fclose(fid);         
        
        fileName = ['Output/drex_',deformationSymmetry,'_info.txt'];
        fid = fopen(fileName,'w');
        fprintf(fid,'DREX Olivine Simulation, Version %3.2f\n',version);
        fprintf(fid,'Model and Analytic Finite-strain Solutions\n');
        fprintf(fid,'**************************\n');
        fprintf(fid,'no deformation');
        fclose(fid);
         
        return
               
    otherwise
        error('specify correct deformationSymmetry');

end 


% check that F has no volume strain
if abs(det(F) - 1) > eps
    error('DREX.m: F must be defined as isochoric (no volume strain)')
end       

%% Decompose Deformation Gradient tensor into flow parameters

Flow.L  = logm(F); % F = expm(Lt);   
Flow.dt = 1/Flow.tSteps; % time is set to 1.

     Flow.D = 0.5*(Flow.L + Flow.L'); % stretching tensor, or rate of deformation tensor
      evals = eig(Flow.D);
Flow.epsnot = max(abs(evals(:)));
Flow.tStar  = Flow.epsnot*Flow.dt*Flow.tSteps;


%% Initialize LPO
%... inital deformation gradient tensor
LPO.deformationGradient = eye(3); 

%... levi-civita epsilon_ijk
Grain.alt = zeros(3,3,3);
Grain.alt(1,2,3) =  1;
Grain.alt(2,3,1) =  1;
Grain.alt(3,1,2) =  1;
Grain.alt(1,3,2) = -1;
Grain.alt(2,1,3) = -1;
Grain.alt(3,2,1) = -1;

% generate initially random texture
LPO.Init.eulerAngles = generaterandomLPO(Grain);

%... assign equal vol fraction to each grain
LPO.Init.volumeFraction   = (1/Grain.nGrains)*ones(Grain.nGrains,1); 
LPO.Init.directionCosines = euler2orientationmatrix(LPO.Init.eulerAngles); 

%% FSE and LPO numerical calculation

% cacluate finite strain, odf, and direction cosines 
[LPO] = strain(Grain,Flow,LPO);


% check if nan
if any(isnan(LPO.Final.volumeFraction));
    warning('At least one ODF Weight is NaN');
    keyboard
end

%... convert direction cosines back to euler angles
LPO.Final.eulerAngles = orientationmatrix2euler(LPO.Final.directionCosines); 


%% Write euler angles to file

writelpotofile(LPO,Flow,Grain);

%% Decompose numerical and analytical deformation gradient tensors

% Decompose analytical (input) deformation gradient tensor
%... left Cauchy Green Tensor
B = F*F';
[~,eval] = eig(B);
eval = sqrt(diag(eval));

%... calculate strain ratios 
[Sx,idx(1)] = max(eval(:));
[Sz,idx(2)] = min(eval(:));
idx(3) = setdiff(1:3,idx);
Sy = eval(idx(3));
clear idx

ShearXY = B(1,2);
% ShearXZ = B(1,3);
% ShearYZ = B(2,3);

% Decompose numerical deformation gradient tensor
fse = LPO.Final.deformationGradient;

%... left cauchy green tensor
Bij = fse*fse';

%... simple shear magnitude = tan(Gam) = LSij(1,3) for simple shear velocity
ShearXYNumerical = Bij(1,2);

%... calculate strain ratios 
[~,eval] = eig(Bij);
eval = sqrt(diag(eval)); 
[SxNumerical,idx(1)] = max(eval);
[SzNumerical,idx(2)] = min(eval);
idx(3) = setdiff(1:3,idx);
SyNumerical = eval(idx(3));

%% Print Output information to file and screen
% set global error flag
global rotationMatrixWarningFlag detTolerance

fileName = ['Output/drex_',Flow.deformationSymmetry,'_info.txt'];
fid = fopen(fileName,'w');
fid = [1,fid];
for i = 1:2
    iFid = fid(i);
    fprintf(iFid,'D-rex Olivine Simulation, Version %s\n',version);
    fprintf(iFid,'Model and Analytic Finite-strain Solutions\n');
    fprintf(iFid,'**************************\n');
    fprintf(iFid,'    Numerical Solution \t Analytical Solution\n');
    fprintf(iFid, 'Shear Strain: \t%f %f \n', ShearXYNumerical, ShearXY);
    fprintf(iFid, 'Sx: %f %f\n',SxNumerical, Sx);
    fprintf(iFid, 'Sy: %f %f\n',SyNumerical, Sy);
    fprintf(iFid, 'Sz: %f %f\n',SzNumerical, Sz);
    fprintf(iFid, 'F =\n %f %f %f\n %f %f %f \n %f %f %f\n',LPO.Final.deformationGradient);
    fprintf(iFid, '******** Parameters **********\n');
    fprintf(iFid, 'esp_0:\t %f\n',   Flow.epsnot);
    fprintf(iFid, 'tSteps:\t %f\n',  Flow.tSteps);
    fprintf(iFid, 't*:\t %f\n',      Flow.tStar);
    fprintf(iFid, 'dt:\t %f\n',      Flow.dt);
    fprintf(iFid, 'nGrains: %i\n',            Grain.nGrains);
    fprintf(iFid, 'pctOlivine: %f \n',        Grain.pctOli);
    fprintf(iFid, 'tau:\t %f %f %f %3.2e \n', Grain.tau);
    fprintf(iFid, 'mob:\t %f \n',             Grain.mob);
    fprintf(iFid, 'chi:\t %f \n',             Grain.chi);
    fprintf(iFid, 'lambda:\t %f \n',          Grain.lambda);
    fprintf(iFid, 'stressExp: %f \n',         Grain.stressExp);
    if rotationMatrixWarningFlag
        % print warning message
        fprintf(iFid, ['at least 1 rotation matrix or orientation matrix is not orthogonal to within tolerance of %f\n',...
                'nearest orthonormal matrix has been substituted for these cases\n******\n\n\n'],detTolerance);
    end

end
fclose(fid(2));
%% Plot using contourstereogram

plotpolefigures(Flow)

%%
% keyboard
end

%% Subfunctions: Output

function [] = writelpotofile(LPO,Flow,Grain)
% Write euler angles to file. All files are written to the Output/
% subfolder. Currently set to write the following files:
% 1. _volfraction.txt: euler angles and volume fraction
% 2. _unweighted.txt:  all non-zero euler angles (one point per grain)
% 3. _volweghted.txt:  euler angles weighted by volume fraction
% 4. __inversevolweighted.txt: euler angles weighted by inverse volume fraction

%% Check inputs
narginchk(3,3);
validateattributes(LPO,{'struct'},{});
validateattributes(Flow,{'struct'},{});
validateattributes(Grain,{'struct'},{});

%% Begin subfunction

% 1. write euler angles and volume fraction to file
fileName = ['Output/drex_',Flow.deformationSymmetry,'_volfraction.txt'];
fid = fopen(fileName,'w');
for i = 1:Grain.nGrains
    fprintf(fid,'%7.3f %7.3f %7.3f %01.5f\n', LPO.Final.eulerAngles(i,:)*180/pi, LPO.Final.volumeFraction(i));
end
fclose(fid); 

% 2. Write all non-zero euler angles to file
fileName = ['Output/drex_',Flow.deformationSymmetry,'_unweighted.txt'];
fid = fopen(fileName,'w');
for i = 1:Grain.nGrains
    if LPO.Final.volumeFraction(i) > 0
        fprintf(fid,'%f %f %f\n', LPO.Final.eulerAngles(i,:)*180/pi);
    end
end
fclose(fid); 

% 3. Write volume weighted euler angles, using random draws to convert odf
% to a discrete number of orientations, weighted by volume

%... sort by vol fraction
[volFracSorted, idxVolFrac] = sort(LPO.Final.volumeFraction);
eulerAnglesSorted = LPO.Final.eulerAngles(idxVolFrac,:);

%... generate cumulative weight
cumWeight = cumsum(volFracSorted);

%... generate random indices
q = qrandstream('halton', 3, 'Skip',1e3, 'Leap',1e2);
idxGrain = qrand(q,Grain.nGrains);

%... find correct eulerAngle and write to file
fileName = ['Output/drex_',Flow.deformationSymmetry,'_volweighted.txt'];
fid = fopen(fileName,'w');
for i = 1:Grain.nGrains
    % find the maximum cumWeight that is less than the random value.
    %... the euler angle index is +1. For example, if the idxGrain(i) < cumWeight(1), the index should be 1 not zero) 
    cumWeightIdx = numel(cumWeight(cumWeight <= idxGrain(i)))+1;
    fprintf(fid,'%f %f %f\n', eulerAnglesSorted(cumWeightIdx,:)*180/pi);
end
fclose(fid);


% 4. Convert odf to a discrete number of orientations, weighted inversely
% by volume, not including grains that have zero volume.

%... sort by vol fraction
volFracNonZeroIdx = find(LPO.Final.volumeFraction > 0);
    volFracNonZero = LPO.Final.volumeFraction(volFracNonZeroIdx);
eulerAnglesNonZero = LPO.Final.eulerAngles(volFracNonZeroIdx);

[volFracSorted, idxVolFrac] = sort(1-volFracNonZero);
    volFracSorted = volFracSorted/sum(volFracSorted);
eulerAnglesSorted = eulerAnglesNonZero(idxVolFrac,:);

%... generate cumulative weight
cumWeight = cumsum(volFracSorted);

% generate random indices
idxGrain = rand(1,Grain.nGrains);

%... find correct eulerAngle and write to fil
fileName = ['Output/drex_',Flow.deformationSymmetry,'_inversevolweighted.txt'];
fid = fopen(fileName,'w');
for i = 1:Grain.nGrains
    % find the maximum cumWeight that is less than the random value.
    %... the euler angle index is +1. For example, if the idxGrain(i) < cumWeight(1), the index should be 1 not zero) 
    cumWeightIdx = numel(cumWeight(cumWeight <= idxGrain(i)))+1;
    fprintf(fid,'%f %f %f\n', eulerAnglesSorted(cumWeightIdx,:)*180/pi);
end
fclose(fid);



end

function [] = plotpolefigures(Flow)
% This subfunction plots the volume weighted euler angles and saves the
% figure using export fig if the export_fig function is available. If the
% export_fig function is not availalbe, the code does not save the figure.
% The volume weighted euler angles are the only plot included by default.
% Code appropriate for plotting the additional outputs is included, but
% commented out.

%% Check inputs
narginchk(1,1);
validateattributes(Flow,{'struct'},{});
validateattributes(Flow.deformationSymmetry,{'char'},{});

%% Begin subfunction

% volume weighted
fileName = ['Output/drex_',Flow.deformationSymmetry,'_volweighted.txt'];
fid = fopen(fileName);
angles = textscan(fid,'%f %f %f %*[^\r\n]','HeaderLines',0);
eulerAngles = [angles{1},angles{2},angles{3}]*pi/180;
clear angles;   
fclose(fid);  

close all
figure(2); clf
[hFig] = contourpolefigures(eulerAngles,'olivine','Gaussian',2);
fileName = ['Output/drex_',Flow.deformationSymmetry,'_volweighted'];

if exist('export_fig','file') == 2
    % save figure
    export_fig(fileName,'-png');
else
    % do nothing
end

% % unweighted
% fileName = ['Output/drex_',Flow.deformationSymmetry,'_unweighted.txt'];
% fid = fopen(fileName);
% angles = textscan(fid,'%f %f %f %*[^\r\n]','HeaderLines',0);
% eulerAngles = [angles{1},angles{2},angles{3}]*pi/180;
% clear angles;   
% fclose(fid);  
% 
% close all
% figure(1); clf
% [hFig] = contourpolefigures(eulerAngles,'olivine','Gaussian',1);
% fileName = ['Output/drex_',Flow.deformationSymmetry,'_unweighted'];
% export_fig(fileName,'-png');


% % inverse volume weighted
% fileName = ['Output/drex_',Flow.deformationSymmetry,'_inversevolweighted.txt'];
% fid = fopen(fileName);
% angles = textscan(fid,'%f %f %f %*[^\r\n]','HeaderLines',0);
% eulerAngles = [angles{1},angles{2},angles{3}]*pi/180;
% clear angles;   
% fclose(fid);  
% 
% close all
% figure(1); clf
% [hFig] = contourpolefigures(eulerAngles,'olivine','Gaussian',1);
% export_fig(fileName,'-png');

end

%% Subfunctions: Calculations

function [eulerAngles] = generaterandomLPO(Grain)
% this function generates a uniform random LPO for use as an initial
% texture in DREX. output lpo is a size3^3 x 3 x 3 matrix. Each index in
% the first dimension is a grain, the second two dimensions give the
% orientation tensor. This uses the method given on p 90 in Morawiec, A.
% (2003). Orientations and rotations. Springer-Verlag.


%% Check inputs
narginchk(1,1);
validateattributes(Grain,{'struct'},{})
validateattributes(Grain.nGrains,{'numeric'},{'scalar'});
%% Begin Function


% updated to use quasi random stream.
nGrains = Grain.nGrains;
q = qrandstream('halton', 3, 'Skip',1e4, 'Leap',1e3);

%... get random values
randVals = qrand(q,nGrains);
eulerAngles(:,1) = 2*pi*randVals(:,1);        
eulerAngles(:,3) = 2*pi*randVals(:,2);       
eulerAngles(:,2) = acos(2*randVals(:,3)-1);   



end

function [odfi] = normalizedvolumefractions(odfi)
% this subfunction checks the volume fraction for impossible values


%% Check inputs
narginchk(1,1);
validateattributes(odfi,{'numeric'},{'ncols',1});


%% Begin Function
% odf weights cannot be negative
odfi(odfi<0) = 0;

% normalize to unit volume
odfi = odfi/sum(odfi(:));

end

function [acsi,varargout] = checkdirectioncosines(acsi)
% Check that all acsi are rotation matrices, to within a set tolerance.
% Includes an output flag later used to print a warning message
%% check inputs
validateattributes(acsi,{'numeric'},{'ncols',9});

% set global error flag
global rotationMatrixWarningFlag detTolerance

%% Begin function
I = eye(3);
nGrains = size(acsi,1);
detTolerance = 1e-6;

for iGrain = 1:nGrains
    a = reshape(acsi(iGrain,:),3,3);
    detA = det(a);
%     fprintf('%e\n',detA);
    if detA-1 > detTolerance
        % abs value is not needed here, because det(rotation) = +1 only.
        % assign warning to print
        rotationMatrixWarningFlag = 1;

        % find nearest orthonormal matrix using singular value method
        [U,~,V] = svd(a);
        a = U*I*V';        
        acsi(iGrain,:) = a(:);
    end
    
end

varargout{1} = rotationMatrixWarningFlag;



end

function [LPO] = strain(Grain,Flow,LPO)
% This function calculates the LPO evolution using a 4th order RK routine.
% The function is setup only for a constant velocity gradient tensor.
% 
% Modified from DREX fortran subroutine STRAIN.


%% Check input arguments
narginchk(3,3);
validateattributes(Grain,{'struct'},{});
validateattributes(Flow, {'struct'},{});
validateattributes(LPO,  {'struct'},{});

%% Begin function

% Parse initial values
deformationGradient = LPO.deformationGradient;
volumeFraction      = LPO.Init.volumeFraction;
directionCosines    = LPO.Init.directionCosines;

% parse Flow parameters
     L = Flow.L; 
epsnot = Flow.epsnot;

    dt = Flow.dt;
nSteps = Flow.tSteps;

% assign initial values for LPO
iDeformationGradient = deformationGradient;
iVolumeFraction      = volumeFraction;
iDirectionCosines    = directionCosines;

dispstat(sprintf('Looping over timesteps:\t'),'keepthis','timestamp');
for iStep = 1:nSteps % loop over all time increments
    dispstat(sprintf('Progress: %03.2f',100*iStep/nSteps),'timestamp');
    
    %% RK Step 1
    [dotacs,dotodf] = calculatederivatives(iDirectionCosines,iVolumeFraction,Grain,Flow);
    
    
    kfse1 = L*iDeformationGradient*dt;
    kodf1 = epsnot*dotodf*dt; % epsnot*dt = dtstar
     kac1 = epsnot*dotacs*dt;
    
    %% RK Step 2
    %... f(xn+0.5dt, yn+0.5*k1)
    iDeformationGradient = deformationGradient + 0.5*kfse1;
    iVolumeFraction      = volumeFraction      + 0.5*kodf1;
    iDirectionCosines    = directionCosines    + 0.5* kac1;

    % check direction cosines and volume fractions
%     iVolumeFraction = normalizedvolumefractions(iVolumeFraction);
    
    % calculate derivatives (evaluate f())    
    [dotacs,dotodf] = calculatederivatives(iDirectionCosines,iVolumeFraction,Grain,Flow);

    % calculate intermediate for RK
    kfse2 = L*iDeformationGradient*dt;
    kodf2 = epsnot*dotodf*dt; % epsnot*dt = dtstar
     kac2 = epsnot*dotacs*dt;


    %% RK Step 3
    %... f(xn+0.5h, yn_0.5k2)
    iDeformationGradient = deformationGradient + 0.5*kfse2;
    iVolumeFraction      = volumeFraction      + 0.5*kodf2;
    iDirectionCosines    = directionCosines    + 0.5* kac2;
    
    % check direction cosines and volume fractions
%     iVolumeFraction = normalizedvolumefractions(iVolumeFraction);
%     [iDirectionCosines,iVolumeFraction] = checkvolumesandcosines(iDirectionCosines,iVolumeFraction);
  
    % calculate derivatives (evaluate f())
    [dotacs,dotodf] = calculatederivatives(iDirectionCosines,iVolumeFraction,Grain,Flow);
      
    % calculate intermediate for RK
    kfse3 = L*iDeformationGradient*dt;
    kodf3 = epsnot*dotodf*dt; 
     kac3 = epsnot*dotacs*dt;
     
    
    %% RK: Step 4
    %... f(xn+h, yn+k3)
    iDeformationGradient = deformationGradient + kfse3;
    iVolumeFraction      = volumeFraction      + kodf3;
    iDirectionCosines    = directionCosines    +  kac3;
    
    % check direction cosines and volume fractions
%     iVolumeFraction = normalizedvolumefractions(iVolumeFraction);
%     [iDirectionCosines,iVolumeFraction] = checkvolumesandcosines(iDirectionCosines,iVolumeFraction);

    % calculate derivatives (evaluate f())
    [dotacs,dotodf] = calculatederivatives(iDirectionCosines,iVolumeFraction,Grain,Flow);
    
    % calculate intermediate for RK
    kfse4 = L*iDeformationGradient*dt;
    kodf4 = epsnot*dotodf*dt;
     kac4 = epsnot*dotacs*dt;
    
    %% RK: Combine all intermediate steps


    deformationGradient = deformationGradient + (0.5*kfse1 + kfse2 + kfse3 + 0.5*kfse4)/3;
    volumeFraction      = volumeFraction      + (0.5*kodf1 + kodf2 + kodf3 + 0.5*kodf4)/3;
                               iGrainRotation = (0.5* kac1 +  kac2 +  kac3 + 0.5* kac4)/3;
    
    % check and update direction cosines:
    iGrainRotation    = checkdirectioncosines(iGrainRotation);      
    directionCosines  = directionCosines + iGrainRotation;    
    directionCosines  = checkdirectioncosines(directionCosines);      
    
    
    %... check for unrealistic direction cosines and volume fractions
    volumeFraction = normalizedvolumefractions(volumeFraction);
    
    % check for nans in direction cosines
    if any(isnan(iDirectionCosines(:)))
        error('strain.m: direction cosines contain NaN values')
    end
    
    % Update tensors for next iteration
    iDeformationGradient = deformationGradient;
    iVolumeFraction      = volumeFraction;
    iDirectionCosines    = directionCosines;    
    

end
dispstat('finished','timestamp');


% parse output
LPO.Final.deformationGradient   = deformationGradient;
LPO.Final.volumeFraction        = volumeFraction;
LPO.Final.directionCosines      = directionCosines;


end

function [dotacs,dotodf] = calculatederivatives(acsi,odfi,Grain,Flow)
% This function caclulates the derivative of the direction cosine matrix
% (dotacs) and the derivative of the odf (dotodf) for each crystal
% orientation. These derivatives are used in the RK integration scheme (see
% strain.m).


%% Check input arguments
narginchk(4,4);

validateattributes(acsi, {'numeric'},{'ncols',9});
validateattributes(odfi, {'numeric'},{'ncols',1});
validateattributes(Grain,{'struct'} ,{});
validateattributes(Flow, {'struct'} ,{});

%% Begin function

% matlab version of subroutine deriv (DREX line 1251)

% Parse inputs 
nGrains    = Grain.nGrains;
L          = Flow.L;
D          = Flow.D;
epsnot     = Flow.epsnot;
tau        = Grain.tau;
stressExp  = Grain.stressExp;
alt        = Grain.alt;
chi        = Grain.chi;
pctOlivine = Grain.pctOli;
Mob        = Grain.mob;
lambda     = Grain.lambda;

% keyboard
% non-dimensionalize velocity gradient and streching tensors
L = L/epsnot;
D = D/epsnot;

% initialize variables
dotacs = nan(nGrains,9);
rt = nan(nGrains,1);
for iGrain = 1:nGrains
    
    % calculate invariants e_{pr} T_{pr} for the 4 slip systems of Ol
    % ... note that this is the function of the current orientation, acsi
    
    % reshape acsi back to 3 x 3
    iGrainAcsi = [acsi(iGrain,1), acsi(iGrain,2), acsi(iGrain,3);
                  acsi(iGrain,4), acsi(iGrain,5), acsi(iGrain,6);
                  acsi(iGrain,7), acsi(iGrain,8), acsi(iGrain,9)];
    
    % I = li * nj * Dij (Kaminski and Ribe, 2001)
    bigi = [0,0,0,0];
    for i1 = 1:3
        for i2 = 1:3
            bigi(1) = bigi(1) + D(i1,i2)*iGrainAcsi(1,i1)*iGrainAcsi(2,i2);
            bigi(2) = bigi(2) + D(i1,i2)*iGrainAcsi(1,i1)*iGrainAcsi(3,i2);
            bigi(3) = bigi(3) + D(i1,i2)*iGrainAcsi(3,i1)*iGrainAcsi(2,i2);
            bigi(4) = bigi(4) + D(i1,i2)*iGrainAcsi(3,i1)*iGrainAcsi(1,i2);
        end
    end
    
    % quotient, I/tau
    q = bigi./tau;
    
    % Reorder quotients I/tau according to absolute magnitude
    qab = abs(q); 
    [~,ti] = sort(qab,'descend');
    imax = ti(1);
    iint = ti(2);
    imin = ti(3);
    iinac= ti(4);  

    % calculate weighting factors gam_s relative to value gam_i for which I/tau is largest 
    %... gam is the ratio of the strain between softest slip system and slip
    % system s
    gam = nan(4,1);
    gam(imax) = 1;
    
    rat  = tau(imax)/bigi(imax);
    qint = rat*bigi(iint)/tau(iint);
    qmin = rat*bigi(imin)/tau(imin);
    sn1  = stressExp-1;
    
    gam(iint) = qint*(abs(qint))^sn1; % Eq 5
    gam(imin) = qmin*(abs(qmin))^sn1;
    gam(iinac)= 0; % inactive slip system
    
    % calculate g, the slip tensor
    g = nan(3,3);
    for i1 = 1:3
        for i2 = 1:3
            g(i1,i2) = 2*(gam(1)*iGrainAcsi(1,i1)*iGrainAcsi(2,i2) + ...
						  gam(2)*iGrainAcsi(1,i1)*iGrainAcsi(3,i2) + ...
						  gam(3)*iGrainAcsi(3,i1)*iGrainAcsi(2,i2) + ...
						  gam(4)*iGrainAcsi(3,i1)*iGrainAcsi(1,i2));
        end
    end
    
    % calculate strain rate of softest slip system
    R1 = 0; R2 = 0;
    for j = 1:3
        i2 = j+2;
        if (i2 > 3)
            i2 = i2-3;
        end
        
        R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j));
		R2 = R2 - (g(j,i2)-g(i2,j))*(L(j,i2)-L(i2,j));
    
        for k = 1:3
             R1 = R1 + 2*g(j,k)*g(j,k);
			 R2 = R2 + 2*L(j,k)*g(j,k);
        end
    end
    
    gam0 = R2/R1;
    
    
    % dislocation density calculation (correct from fortran version)
    rt1=tau(imax)^(1.5-stressExp)*abs(gam(imax)*gam0)^(1.5/stressExp);
    rt2=tau(iint)^(1.5-stressExp)*abs(gam(iint)*gam0)^(1.5/stressExp);
    rt3=tau(imin)^(1.5-stressExp)*abs(gam(imin)*gam0)^(1.5/stressExp);
    rt4=tau(iinac)^(1.5-stressExp)*abs(gam(iinac)*gam0)^(1.5/stressExp);
    
    rt(iGrain) = rt1*exp(-lambda*rt1^2) +... 
                 rt2*exp(-lambda*rt2^2) +... 
                 rt3*exp(-lambda*rt3^2) +...
                 rt4*exp(-lambda*rt4^2);   
 

    % calculation of rotation rate (this is the spin vector, rot = del
    % cross v - the spin vector for the crystal)
    rot(3) = (L(2,1)-L(1,2))/2-((g(2,1)-g(1,2))/2)*gam0;
    rot(2) = (L(1,3)-L(3,1))/2-((g(1,3)-g(3,1))/2)*gam0;
    rot(1) = (L(3,2)-L(2,3))/2-((g(3,2)-g(2,3))/2)*gam0;
    
    % derivative of the dir cosine matrix
    iGrainDotAcsi = zeros(3,3);
    for i1 = 1:3
        for i2 = 1:3
            for i3 = 1:3
                for i4 = 1:3
                    iGrainDotAcsi(i1,i2) = iGrainDotAcsi(i1,i2)+alt(i2,i3,i4)*iGrainAcsi(i1,i4)*rot(i3);
                end
            end
        end
    end
    
    % grain boundary sliding for small grains
    if (odfi(iGrain) < chi/nGrains)
        % sliding grains do not rotate or accumulate dislocation energy
        iGrainDotAcsi =  zeros(3,3);
        rt(iGrain) = 0;
    end
    iGrainDotAcsi = iGrainDotAcsi';
    dotacs(iGrain,1:9) = iGrainDotAcsi(:);
    
end % iGrain

% volume averaged energy
Emean = sum(odfi(:).*rt(:)); % sometimes bigi = zeros
    
% change vol fraction by grain boundary migration and energy
dotodf = pctOlivine * Mob * odfi(:) .* (Emean-rt(:));

% make col vec   
dotodf = dotodf(:);  




end


