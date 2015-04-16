function CrystalDirections = kambcontour(CrystalDirections,LambertProj)


% parse inputs
gridVectors = [LambertProj.x(:), LambertProj.y(:), LambertProj.z(:)];
[nSpherePoints,~] = size(LambertProj.X);


% Calculate size of counting circle (window parameter)
nData = CrystalDirections.nData;
m     = CrystalDirections.multiplicity;
theta = acos(1-9/(m*(nData+9)));

%... for extremely large numbers of directions (say, over 150,000) this
%vectorized version can overload ram and become extremely slow. Include a
%contingency plan
if nData*m < 150000
    % vectorized version (much faster)
    data = CrystalDirections.unitVectors;
    
    % Calculate the angle between the spherePoint and each of the data points
    %... vectorized version (doesn't use cross-product)
    cosAng = gridVectors*data';
    clear gridVectors    
    
    fprintf(1,'Vectorized angular distance calculation for %s \n',CrystalDirections.name);
    angDist = acos(abs(cosAng)); % absolute value ensures acute angle is always returned
    clear cosAng
    
    % count distances less than theta
    counts = sum(angDist<=theta,2); clear angDist
    
else
    % large set of angles (e.g. m poles in quartz has 6xdata). Use loops as
    % matrices are likely much too large for RAM

    counts = zeros(size(gridVectors(:,1)));
    dispstat(sprintf('Calculating angular distances for %s: \t',CrystalDirections.name),'keepprev','keepthis','timestamp');
    for i = 1:nSpherePoints
        dispstat(100*i/nSpherePoints);

       % operate on a subset of sphere points
       idx = (i-1)*nSpherePoints + 1:(i-1)*nSpherePoints + nSpherePoints;
       
       % calculate dot product for this subset
        gridV = gridVectors(idx,:); % vector of grid point
       cosAng = gridV*CrystalDirections.unitVectors'; % dot product
       cosAng = acos(abs(cosAng)); % get angular distance
       
       counts(idx) = sum(cosAng <= theta,2); % number inside counting circle
       
    end
    dispstat('done');

    clear cosAng

end

% reshape counts back to square matrix
counts = reshape(counts,nSpherePoints,nSpherePoints);


% normalize counts by sigma (std deviation), assuming uniformly distributed
% binomial trials
areaFrac = 9/(nData+9);
   sigma = sqrt(nData*areaFrac*(1-areaFrac));
  counts = counts/(3*sigma);


% Parse output

CrystalDirections.parameter = theta;
CrystalDirections.paramName = sprintf('Window = %02.2f%c',theta,char(176));
CrystalDirections.counts    = counts;
CrystalDirections.minCounts = min(counts(:));
CrystalDirections.maxCounts = max(counts(:));
