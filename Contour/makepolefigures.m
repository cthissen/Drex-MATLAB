function [hFig] = makepolefigures(CrystalDirections,mineral,SphereProj,hFig,ColorOpts)

% Parse Values common to each axis (coordinate frame)
X = SphereProj.X;
Y = SphereProj.Y;
R = SphereProj.R; % radius of circle in X and Y

% Setup Figure Properties
hFig.Units    = 'centimeters';
hFig.Position = [0,0,19,10]; % 2 column figure width is 19 cm 
hFig.Name     = 'Pole Figures';
movegui(hFig,'northwest');
fontSize = 18;

% Calculate pole figure spacing for three pole figures + 1 infobox
figBox = 19;
nAxesHz = 3; % 4 horizontal subaxes
subaxesSpacing = 0.25*(1/figBox); % 1/4 cm spacing between axes
totalSpacing = (nAxesHz+1)*subaxesSpacing; % total spacing (normalized)
subaxesWidth = (1-totalSpacing)/nAxesHz;

colorRamp0     = ColorOpts.colorRamp0;
nTics          = ColorOpts.nTics;
zeroValue      = ColorOpts.centerVal;
factor         = ColorOpts.factor;

% Now make pole figure plot
switch mineral
    case 'olivine'

    %... create common colormap    
    minCounts = min([CrystalDirections(:).minCounts]);
                     
    [maxCounts,idx] = max([CrystalDirections(:).maxCounts]);
                 

     
     %... modify counts for proper scaling
     mask = zeros(size(X));
     mask(sqrt(X.^2 + Y.^2) >= R +0.025) = NaN; % mask values outside net
     mask(sqrt(X.^2 + Y.^2) <  R +0.025) = 1;      
     counts = CrystalDirections(idx).counts.*mask;
     counts = counts(~isnan(counts(:)));
%     [colorRamp1,tics] = cmapscale(log(counts+1e-4),colorRamp0,factor,zeroValue,nTics); 
    [colorRamp1,tics] = cmapscale(counts,colorRamp0,factor,zeroValue,nTics); 


    % Plot Stereograms (1 row only)
    topRow = 0;
    for idx = 1:3        
        hSub(idx) = axes;
        leftPos = idx*subaxesSpacing + (idx-1)*subaxesWidth;
        hSub(idx).Position = [leftPos,topRow,subaxesWidth,1];

        if exist('padarray','builtin')
            % use gcolor for slightly betterplotting
            h = gcolor(X,Y,CrystalDirections(idx).counts.*mask); shading flat
        else
            % use pcolor to plot
            h=pcolor(X,Y,mask.*CrystalDirections(idx).counts); shading interp
        end
        hold on
        contour(X,Y,CrystalDirections(idx).counts.*mask,[2,2],'Color',[0.5,0.5,0.5],'LineWidth',2);

        set(gca,'YDir','normal', 'Color', 'none');
        set(h,'alphadata',~isnan(CrystalDirections(idx).counts.*mask));
            colormap(colorRamp1);
            hold on
            polenet(R);
            axis equal off
%             hTitle = title(CrystalDirections(idx).name,'FontSize',fontSize);   
%             hTitle.Position = [0,1.5,0];
        hSub(idx).CLim = [minCounts,maxCounts];
    end

    % Plot the InfoBox
    %... create invisible axes the full size of the figure (makes locating textbox easier)
    hAxis = axes('units','normalized','pos',[0 0 1 1],...
                 'visible','off','handlevisibility','off');
             
    % Colorbar
    leftPos = 3*subaxesSpacing + 2*subaxesWidth + 0.6*subaxesWidth;
    colormap(hAxis,colorRamp0)
    hCbar = colorbar(hAxis);

    cbarWidth = 0.25*subaxesWidth;
    hCbar.Location = 'northoutside'; % makes colorbar hz
    hCbar.Position = [leftPos,0.15,cbarWidth,0.03];
    hCbar.FontSize = 14;
    hCbar.Title.String = 'MUD';
    hCbar.Title.Units = 'normalized';
    hCbar.Title.Position = [0.5,-4,0];    

    set(hAxis,'CLim',[0,1])    
    hCbar.Ticks = tics(:,1);
    hCbar.TickLabels = sprintf('%3.0f\n',tics(:,2));
        
    
    % Textbox
    leftPos = 3*subaxesSpacing + 2*subaxesWidth;
    nData = CrystalDirections(1).nData;
    theta = CrystalDirections(1).parameter;
  
    
    stringLimited = sprintf(['N = %i \nMax: %03.2f'],...
        nData,maxCounts);
    
    text(leftPos,0.15,stringLimited,'parent',hAxis,'FontSize',13)

    
    case 'quartz'
    hFig.Position = [0,0,19,16]; % 2 column figure width is 19 cm 

    
    
    %... create common colormap    
    minCounts = min([CrystalDirections(:).minCounts]);
                     
    [maxCounts,idx] = max([CrystalDirections(:).maxCounts]);
                 
     %... modify counts for proper scaling
     mask = zeros(size(X));
     mask(sqrt(X.^2 + Y.^2) > sqrt(2)+0.025) = NaN; % mask values outside net
     mask(sqrt(X.^2 + Y.^2) <  sqrt(2)+0.025) = 1;      
     counts = CrystalDirections(idx).counts.*mask;
     counts = counts(~isnan(counts(:)));
    [colorRamp1,tics] = cmapscale(counts,colorRamp0,factor,zeroValue,nTics); 

    
    
    % Plot Stereograms (TOP ROW)
    topRow = 0.25;
    for idx = 1:3        
        hSub(idx) = axes;
        leftPos = idx*subaxesSpacing + (idx-1)*subaxesWidth;
        hSub(idx).Position = [leftPos,topRow,subaxesWidth,1];
        pcolor(X,Y,mask.*CrystalDirections(idx).counts); shading interp
            colormap(colorRamp1);
            hold on
            contour(X,Y,mask.*CrystalDirections(idx).counts,[2,2],'Color',[0.5,0.5,0.5],'LineWidth',2);            
            polenet(sqrt(2));
            axis equal off
            hTitle = title(CrystalDirections(idx).name,'FontSize',fontSize);   
            hTitle.Position = [0,1.5,0];
        hSub(idx).CLim = [minCounts,maxCounts];
    end
    
    % Plot Stereograms (Bottom ROW)
    for idx = 4:5
    bottomRow = -0.25;
        hSub(idx) = axes;
        leftPos = (idx-3)*subaxesSpacing + (idx-4)*subaxesWidth + 0*subaxesWidth;
        hSub(idx).Position = [leftPos,bottomRow,subaxesWidth,1];
%         contourf(X,Y,mask.*CrystalDirections(idx).counts,contourLevels);
        pcolor(X,Y,mask.*CrystalDirections(idx).counts); shading interp
            colormap(colorRamp1);
            hold on
            contour(X,Y,mask.*CrystalDirections(idx).counts,[2,2],'Color',[0.5,0.5,0.5],'LineWidth',2);
            polenet(sqrt(2));
            axis equal off
            hTitle = title(CrystalDirections(idx).name,'FontSize',fontSize);   
            hTitle.Position = [0,1.5,0];
        hSub(idx).CLim = [minCounts,maxCounts];   
    end
    % Plot the InfoBox
    %... create invisible axes the full size of the figure (makes locating textbox easier)
    hAxis = axes('units','normalized','pos',[0 0 1 1],...
                 'visible','off','handlevisibility','off');
             
    %... colorbar
    leftPos = 4*subaxesSpacing + 2*subaxesWidth + 0*subaxesWidth;
    colormap(hAxis,colorRamp0)
    hCbar = colorbar(hAxis);

    cbarWidth = 0.4*subaxesWidth;
    hCbar.Location = 'northoutside'; % makes colorbar hz
    hCbar.Position = [leftPos,0.25,cbarWidth,0.03];
    hCbar.FontSize = 14;
    hCbar.Title.String = 'MUD';
    hCbar.Title.Units = 'normalized';
    hCbar.Title.Position = [0.5,-2.5,0];
    
    set(hAxis,'CLim',[0,1])    
%     hCbar.Limits = [0,1];
    hCbar.Ticks = tics(:,1);
    hCbar.TickLabels = sprintf('%3.0f\n',tics(:,2));
    

    
    
    % Textbox
    nData = CrystalDirections(1).nData;
    theta = CrystalDirections(1).parameter;
    thetam3 = CrystalDirections(3).parameter;
    
    stringFull = sprintf(['N = %i \nMethod: Kamb \nRSD=1/3 \nWindow = %04.2f%c \nEqual',...
        ' Area Projection \nLower Hemisphere'],...
        nData,theta*180/pi,char(176));
    
    stringLimited = sprintf([' N = %i \n [c]:  %s',...
        '\n<a>: %s \n Max: %03.2f'],...
        nData,CrystalDirections(1).paramName,CrystalDirections(2).paramName,maxCounts);
    
    stringLimited = sprintf(' N = %i',...
        nData);

    
    text(leftPos,0.35,stringLimited,'parent',hAxis,'FontSize',13)
    

    otherwise
        error('contourpolefigures.m: Specify quartz or olivine for mineral');
    
end

end