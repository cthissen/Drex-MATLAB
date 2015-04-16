function Directions = gaussiancontour(Directions,SphereProj)
% Calculate a Gaussian smoothed kernel of the input Directions over the
% points given in SphereProj. This function uses the absolute value of the
% cosine of the angle to give the acute angle between the grid points and
% the data.
%
% Directions is a structure that contains the Nx3 unitVectors field.
% Directions also contains nData for the number of data and multiplicity,
% which gives number of equivalent directions. 
% SphereProj is also a structure that contains the x, y, and z
% coordinates of the grid on which the Gaussian smoothing is calculated.
%
% For more information, see: Robin, P. Y. F., & Craig Jowett, E. (1986).
% Computerized density contouring and statistical evaluation of orientation
% data using counting circles and continuous weighting functions.
% Tectonophysics, 121(2), 207-223.
% 
% 
% parse inputs
gridVectors = [SphereProj.x(:), SphereProj.y(:), SphereProj.z(:)];
[nSpherePoints,~] = size(SphereProj.X);

% Calculate window parameter, k
nData = Directions.nData;
    m = Directions.multiplicity;
    k = 2*(1+m*nData/9); % Robin and Jowett, Table 3 pdf pg 9
    
%     k = k/2
%     keyboard
 
maxData = 100000;
% keyboard
if nData*m <= maxData && numel(SphereProj.x) < 3000000
    % vectorized version (much faster)
    
    % Calculate the angle between the spherePoint and each of the data points
    cosAng = gridVectors*Directions.unitVectors';
    clear gridVectors    

    fprintf(1,'Vectorized angular distance calculation for %s \n',Directions.name);

    % calculate the continuous weighting function
    counts = exp(k*(abs(cosAng)-1)); % absolute value ensures acute angle is always returned
    clear cosAng
    
    % w is a nSpherePoints x nData matrix. Sum along the rows to calculate
    % the count value for each grid point
    counts = sum(counts,2);
    
    
else
    % loop for large datasets
    counts = zeros(size(gridVectors(:,1)));
    dispstat(sprintf('Calculating angular distances for %s: \t',Directions.name),'keepthis');
    for i = 1:nSpherePoints
        dispstat(sprintf('Progress %02.1f%%',100*i/nSpherePoints)); 
        
        % operate on a subset of sphere points (partially vectorized)
        idx = (i-1)*nSpherePoints + 1:(i-1)*nSpherePoints + nSpherePoints;
        
        % calculate dot product for this subset
        gridV = gridVectors(idx,:); % vector of grid point
        cosAng = gridV*Directions.unitVectors';
        
        % calculate the continuous weighting function
        tmp = exp(k*(abs(cosAng)-1)); % absolute value ensures acute angle is always returned
        clear cosAng
        
        % w is a nSpherePoints x nData matrix. Sum along the rows to calculate
        % the count value for this
        counts(idx) = sum(tmp,2);
        clear tmp
    end
    dispstat('');
    
    
    
    
end

% calculate standard deviation (eq 13b)
stdDev = sqrt(m*nData*(k/2-1)/(k^2));

% normalize so each MUD is 3 sigma from that expected for a uniform
% distribution
counts = counts/(3*stdDev);

% reshape counts back to square matrix
counts = reshape(counts,nSpherePoints,nSpherePoints);
  
% Parse output
Directions.parameter = stdDev;
Directions.paramName = sprintf('\\sigma = %02.2f',stdDev);
Directions.counts    = counts; 
Directions.minCounts = min(counts(:));
Directions.maxCounts = max(counts(:));


end
