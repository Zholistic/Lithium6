function [ handleG ] = plotErrorLines( xArray, yArray, xError, yError, argXY, fignumber, color, markerType )



if(argXY == 0)
    %only yError:
    handleG = figure(fignumber);
    plot(xArray,yArray,'Color',color,'MarkerSize',5,...
        'MarkerFaceColor',color,...
        'Marker',markerType,...
        'LineStyle','--');
    hold on;
    
    for i=1:length(yArray)
        line([xArray(i) xArray(i)],[(yArray(i) - yError(i)) (yArray(i) + yError(i))], ...
            'Color',color);
    end

    grid on;
    
elseif(argXY == 1)
    
    handleG = figure(fignumber);
    plot(xArray,yArray,'Color',color,'MarkerSize',5,...
        'MarkerFaceColor',color,...
        'Marker',markerType,...
        'LineStyle','--');
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




end