function [ handleG ] = plotErrorLinesFig2aMain( xArray, yArray, xError, yError, argXY, fignumber, color, markerType )
%width = 2.6;     % Width in inches
%height = 1.2;    % Height in inches
width = 6;     % Width in inches
height = 2;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 8;      % Fontsize
lw = 0.5;      % LineWidth
msz = 4;       % MarkerSize
%dmn = -pi:0.001:pi;

handleG = figure(fignumber);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
%set(gca, 'FontSize', fsz, 'LineWidth', alw,'YGrid','on','YTick', 4.2:0.4:6);
%set(gca, 'FontSize', fsz, 'LineWidth', alw,'GridLineStyle','-','YGrid','on','XGrid','on','YMinorTick','on','XMinorTick','on');
set(gca, 'FontSize', fsz, 'LineWidth', alw,'YMinorTick','on','XMinorTick','on'); %<- Set properties
%plot(dmn,f(dmn),'b-','LineWidth',lw,'MarkerSize',msz); %<- Specify plot properites
box on;

if(argXY == 0)
    %only yError:

    hold on;
    
    line([832.2 832.2],[16 22],'Color',[0 0 0],'LineWidth', lw,'LineStyle','-.');    
    
    for i=1:length(yArray)
        line([xArray(i) xArray(i)],[(yArray(i) - yError(i)) (yArray(i) + yError(i))], ...
            'Color',color,'LineWidth', lw);
    end
        %'MarkerFaceColor',[0.701960802078247 0.780392169952393 1],...    
    plot(xArray,yArray,'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',[0.701960802078247 0.780392169952393 1],...
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);
        
    xlim([680 990]);
    ylim([16.5 21.5]);
    %legend('f(x)', 'g(x)', 'f(x)=g(x)', 'Location', 'SouthEast');
    %xlabel('Magnetic Field (G)');
    %ylabel('Transverse Width (\mu m)');
    %title('Transverse Width vs Magnetic Field (N_{\sigma} = 12 thousand)');
        
    
elseif(argXY == 1)
    
    %handleG = figure(fignumber);
    plot(xArray,yArray,'Color',color,'MarkerSize',6,...
        'MarkerFaceColor',[0.701960802078247 0.780392169952393 1],...
        'Marker',markerType,...
        'LineStyle','none');
    hold on;
    
    for i=1:length(yArray)
        line([xArray(i) xArray(i)],[(yArray(i) - yError(i)) (yArray(i) + yError(i))], ...
            'Color',color);
    end
    for i=1:length(xArray)
        line([(xArray(i)-xError(i)) (xArray(i)+xError(i))],[(yArray(i)) (yArray(i))], ...
            'Color',color);
    end

    grid on;
end

%set(gca,'XTick',-3:3);
%set(gca,'YTick',0:10);

% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
%set(gcf,'Renderer','OpenGL')

% Save the file as PNG
print('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2a_Main_v11','-depsc2','-r300');
print('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2a_Main_v11','-dpng','-r300');
saveas(handleG,'C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2a_Main_v11','fig');


end