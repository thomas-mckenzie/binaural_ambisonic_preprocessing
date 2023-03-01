function [] = heatmap_plot(az,el,z)
% Plot spectral difference of a large spherical set of points on a 
% rectangular plot. Need input vectors of azimuth, elevation and the
% spectral difference of the points. 



% fig = figure;%('Position', [40, 40, 650, 500]);
x = az;
y = el;

xlin = linspace(min(x),max(x),180*2);
ylin = linspace(min(y),max(y),90*2);
[X,Y] = meshgrid(xlin,ylin);

Z = griddata(x,y,z,X,Y,'cubic');

surf(X,Y,Z,'EdgeColor','none')
axis tight;
xlabel('Azimuth (°)');
ylabel('Elevation (°)');

% shading interp
view ([0 90])
% s.EdgeColor = 'none';

set(gca, 'ActivePositionProperty' , 'position');

originalSize = get(gca, 'Position');
originalOuterBoundary = get(gca, 'outerposition');

set(gca, 'Position', originalSize);
set(gca, 'outerposition', originalOuterBoundary);

colormap(flipud(parula))
set(gca, 'XDir', 'reverse', 'YTick', -75:75:75, 'XTick', -150:75:150); %%%% should that be -180:75:180   ???
xlim([-180 180]);
ylim([-90 90]);

% m = pbaspect('mode');
% m = 'manual';
pbaspect([2 1 1])
% daspect([1 0.6 1]);

box on
set(gcf, 'Color', 'w');




% 
% if saveFigs == 1
% %     savefig(strcat('Figures301/deltaSpecDiff_order_',num2str(AmbiOrder),'_Both_Ears_ambiDFC_IRs_16020_1024tap_quadBias_',num2str(quadratureSFactor),'.fig'));
%     
%     set(gcf, 'Color', 'w');
% %     export_fig(sprintf('Figures301/deltaSpecDiff_order%d_bothEars_ambiDFC_IRs_16020_1024tap_quadBias%d.pdf',AmbiOrder,quadratureSFactor), '-pdf', '-q101', '-opengl');
%     
% %     
% %     colormap(gray);
% %     set(gcf, 'Color', 'w');
% %     export_fig(sprintf('Figures301/deltaSpecDiff_order%d_bothEars_ambiDFC_IRs_16020_1024tapGRAY.pdf',AmbiOrder), '-pdf', '-q101', '-opengl');
% end




end

