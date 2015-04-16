

% test colormaps

nColors = 256;
cMap(1).map = algae(nColors);
cMap(2).map = virus(nColors);
cMap(3).map = coolwarm(nColors);
cMap(4).map = hklcolor(nColors);
cMap(5).map = rdylbu(nColors);

cMap(1).name = 'algae';
cMap(2).name = 'virus';
cMap(3).name = 'coolwarm';
cMap(4).name = 'hklcolor';
cMap(5).name = 'rdylbu';


close all
hFig = figure(1); clf

nCmap = numel(cMap);
fontSize = 14;

% calculate spacing
vertWidth = 0.05;
vertPos = 1/(nCmap + 3*nCmap*vertWidth);

for i = 1:nCmap;
    h(i) = axes('units','normalized','pos',[0 0 1 1],...
            'visible','off','handlevisibility','off');
    hCbar = colorbar(h(i));
    hCbar.Location = 'southoutside'; % makes colorbar hz
    colormap(h(i),cMap(i).map);
    hCbar.Position = [0.1, i*vertPos, 0.8, vertWidth];
    hCbar.Label.String = cMap(i).name;
    
    hCbar.FontSize = fontSize;
    hCbar.Label.FontSize = fontSize;
end






