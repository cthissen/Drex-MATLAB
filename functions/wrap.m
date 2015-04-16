function angle = wrap(angle, minAngle, span)
%wrap Input angle is wrapped to [minAngle,minAngle+span].
% All angles should be given in either degrees or radians.
% For example, minAngle = -pi/2 and span = pi would convert the input angle 
% to the interval -pi/2 to +pi/2.  
%% Calculation
angle = angle - span.*floor((angle - minAngle)./span);
end
