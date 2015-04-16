function polenet(R,varargin)
% plots background for polenet

optArgs = {true}; % default is plot cross
nArgsIn = find(~cellfun(@isempty,varargin));
optArgs(nArgsIn) = varargin(nArgsIn);
[plotCross] = optArgs{:};

theta = (0:pi/180:2*pi);
% R = 1;

x = R*cos(theta);
y = R*sin(theta);
plot(x,y,'k-','LineWidth',1.5)
hold on


switch plotCross
    case true
        plot([0, 0], [-R, R],'k', [-R, R], [0, 0],'k');
        plot([-R, R], [0, 0],'k', [0, 0], [-R, R],'k');
    otherwise
        % don't plot cross
end


end