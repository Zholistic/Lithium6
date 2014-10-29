function [ handleG ] = plotErrorLines( xArray, yArray, xError, yError, argXY, fignumber, color, colorM, markerType, hozFits, slanFits )
%width = 2.6;     % Width in inches
%height = 1.2;    % Height in inches
width = 2;     % Width in inches
height = 2;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 8;      % Fontsize
lw = 0.5;      % LineWidth
msz = 5;       % MarkerSize
%pixelLength = 2.84e-6; %2.84 um topcam
%dmn = -pi:0.001:pi;
dSet = 1;
%mag = 2;

if(fignumber == 1000)
    dSet = 1;
    %mag = 2;
elseif(fignumber == 1001)
    dSet = 2;
    %mag = 7;
elseif(fignumber == 1002)
    dSet = 3;
    %mag = 10;
end

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
    
    for i=1:length(yArray)
        line([xArray(i) xArray(i)],[(yArray(i) - yError(i)) (yArray(i) + yError(i))], ...
            'Color',color,'LineWidth', lw);
    end
        %'MarkerFaceColor',[0.701960802078247 0.780392169952393 1],...    
    plot(xArray,yArray,'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',colorM,...  
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);
    
    if(dSet ~= 1)
    plot(xArray(hozFits),yArray(hozFits),'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',color,...  
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);  
    end
    
    xlim([680 1000]);
    ylim([4.2 5.6]);
    %legend('f(x)', 'g(x)', 'f(x)=g(x)', 'Location', 'SouthEast');
    %xlabel('Magnetic Field (G)');
    %ylabel('Transverse Width (\mu m)');
    %title('Transverse Width vs Magnetic Field (N_{\sigma} = 12 thousand)');
        
    
elseif(argXY == 1)
    
    %handleG = figure(fignumber);
    hold on;
    for i=1:length(yArray)
        line([xArray(i) xArray(i)],[(yArray(i) - yError(i)) (yArray(i) + yError(i))], ...
            'Color',color,'LineWidth', lw);
    end
    for i=1:length(xArray)
        line([(xArray(i)-xError(i)) (xArray(i)+xError(i))],[(yArray(i)) (yArray(i))], ...
            'Color',color,'LineWidth', lw);
    end
    
    plot(xArray,yArray,'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',colorM,...
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);
    
    if(dSet ~= 1)
    plot(xArray(hozFits),yArray(hozFits),'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',color,...  
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);  

    plot(xArray(slanFits),yArray(slanFits),'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',color,...
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);
    
    end

    
    xlim([0 33800]);
    ylim([15 21]);
    %legend('f(x)', 'g(x)', 'f(x)=g(x)', 'Location', 'SouthEast');
    %xlabel('Magnetic Field (G)');
    %ylabel('Transverse Width (\mu m)');
    %title('Transverse Width vs Magnetic Field (N_{\sigma} = 12 thousand)');
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
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2b_Main' num2str(fignumber) '_v13'],'-depsc2','-r300');
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2b_Main' num2str(fignumber) '_v13'],'-dpng','-r300');

figure(20+fignumber);
    plot(xArray,yArray,'Color',color,'MarkerSize',msz,...
        'MarkerFaceColor',colorM,...  
        'Marker',markerType,...
        'LineStyle','none','LineWidth',lw);
    if(dSet == 1)
    legend '832.2 G';
    elseif(dSet == 2)
        legend '865 G';
    elseif(dSet == 3)
        legend '950 G';
    end
% Save the file as PNG
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2b_Main' num2str(fignumber+20) '_Legend_v3'],'-depsc2','-r300');
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2b_Main' num2str(fignumber+20) '_Legend_v3'],'-dpng','-r300');

end