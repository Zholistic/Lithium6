
%Reprocess images datasets?
reprocess = 1; %TODO

TotalDatasets = 3;

directory = cell(TotalDatasets,1);
datestring = cell(TotalDatasets,1);
fileLocList = cell(TotalDatasets,1);
logSigmaY = cell(TotalDatasets,1);
logSigmaX = cell(TotalDatasets,1);
logNROISum = cell(TotalDatasets,1);

directory{1} = 'C:\Data\140318_833G_Crossover_Measurement_10us_Imagepulse_1_I\';
directory{2} = 'C:\Data\140318_880G_Crossover_Measurement_10us_Imagepulse_1_Isat\';
directory{3} = 'C:\Data\140317_997G_Crossover_Measurement_10us_Imagepulse_0_5_Isat\';

datestring{1} = '140318';
datestring{2} = '140318';
datestring{3} = '140317';

bins(1) = 18;
bins(2) = 18;
bins(3) = 18;

if(reprocess)
for i=1:length(directory(:))
    
    [logSigmaY{i},logSigmaX{i},logNROISum{i}] = TopCamN_CrossoverCalc(directory{i},datestring{i},bins(i));
    
end
end

%24 Bins. Select for those above and below 30,000 atoms.

atomnumber = 27000;
crossover = log10(atomnumber);
crossIndex = [];
binStart = 2; %Default =1, for cutting out low atom number data.


for i=1:length(logNROISum(:))
    a = i;
    for j=1:length(logNROISum{a})
        
        if logNROISum{i}(j) > crossover
            crossIndex(i) = j;
            break;
        end
        
    end
end

%Scan for bad values such as NaN or -Inf
for i=1:length(logNROISum(:))
    a = i;
    for j=1:length(logNROISum{a})
        
        %Set them to zero...
        if isnan(logNROISum{i}(j)) || isinf(logNROISum{i}(j))
            logNROISum{i}(j) = 0;
            logSigmaX{i}(j) = 0;
            logSigmaY{i}(j) = 0;
            disp('Found NaN or Inf! Set it to 0.');
        end
        
    end
end

%for i=1:length(logNROISum(:))
%    sectionA(i) = 1:crossIndex(i);
%    sectionB(i) = crossIndex(i):length(logNROISum{i});
%end

%Now we have all the data loaded and the regions of data on which to plot.

%Fit linear plots to the regions:

fg = @(p,x)(p(1).*x + p(2)); %function to fit with
p0 = [0.25 1];
lb = [0.01 0];
ub = [0.42 3];
%p0 = [0.25 12];
%lb = [0.01 5];
%ub = [0.42 30];

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
xs = 1:length(logNROISum{1}); 

coefsAx = []; coefsAy = []; coefsBx = []; coefsBy = [];

for i=1:length(logNROISum(:))
    
    coefsAx(:,i) = lsqcurvefit(fg,p0,logNROISum{i}(binStart:crossIndex(i)),logSigmaX{i}(binStart:crossIndex(i)),lb,ub,curvefitoptions);
    coefsAy(:,i) = lsqcurvefit(fg,p0,logNROISum{i}(binStart:crossIndex(i)-1),logSigmaY{i}(binStart:crossIndex(i)-1),lb,ub,curvefitoptions);
    coefsBx(:,i) = lsqcurvefit(fg,p0,logNROISum{i}(crossIndex(i)+1:length(logNROISum{i})),logSigmaX{i}(crossIndex(i)+1:length(logNROISum{i})),lb,ub,curvefitoptions);
    coefsBy(:,i) = lsqcurvefit(fg,p0,logNROISum{i}(crossIndex(i)+2:length(logNROISum{i})),logSigmaY{i}(crossIndex(i)+2:length(logNROISum{i})),lb,ub,curvefitoptions);   
    
end

fullxsA = 3.3:0.1:4.7; %Full range to plot against
fullxsB = 4:0.1:5;

%The main plot:
figure('Position',[200, 150, 1600, 800]);
title('Red: Fit up to 40k Atoms. Blue: Fit beyond 40k Atoms.');
for i=1:length(logNROISum(:))
    subplot(2,3,i);
    hold on;
    axis([fullxsA(1) fullxsB(end) 1 1.65]);
    plot(logNROISum{i}(binStart:crossIndex(i)),logSigmaY{i}(binStart:crossIndex(i)),'.','Color','r'); %'xk' original symbol
    plot(logNROISum{i}(crossIndex(i)+1:length(logNROISum{i})),logSigmaY{i}(crossIndex(i)+1:length(logNROISum{i})),'.','Color','b'); %'+k' original symbol
    plot(fullxsA,fg(coefsAy(:,i),fullxsA),'Color','r');
    plot(fullxsB,fg(coefsBy(:,i),fullxsB),'Color','b');
    hold off;
    
    if i==1
        title({'833 Gauss',['Red a = ' num2str(coefsAy(1,i),3) '. Blue a = ' num2str(coefsBy(1,i),3) '.']});
    elseif i==2
        title({'880 Guass',['Red a = ' num2str(coefsAy(1,i),3) '. Blue a = ' num2str(coefsBy(1,i),3) '.']});
    else
        title({'979 Gauss',['Red a = ' num2str(coefsAy(1,i),3) '. Blue a = ' num2str(coefsBy(1,i),3) '.']});
    end

end

for i=1:length(logNROISum(:))
    subplot(2,3,i+length(logNROISum(:)));
    hold on;
    axis([fullxsA(1) fullxsB(end) 1 1.65]);
    plot(logNROISum{i}(1:crossIndex(i)),logSigmaX{i}(1:crossIndex(i)),'xk');
    plot(logNROISum{i}(crossIndex(i)+1:length(logNROISum{i})),logSigmaX{i}(crossIndex(i)+1:length(logNROISum{i})),'+k');
    plot(fullxsA,fg(coefsAx(:,i),fullxsA),'Color','r');
    plot(fullxsB,fg(coefsBx(:,i),fullxsB),'Color','b');
    hold off;
    title(['Red a = ' num2str(coefsAx(1,i),3) '. Blue a = ' num2str(coefsBx(1,i),3) '.']);

end

figure(2)
i=3;
hold on;
%axis([fullxsA(1) fullxsB(end) 1 1.65]);
plot(logNROISum{i}(binStart:crossIndex(i)),logSigmaY{i}(binStart:crossIndex(i)),'.','Color','r'); %'xk' original symbol
plot(logNROISum{i}(crossIndex(i)+1:length(logNROISum{i})),logSigmaY{i}(crossIndex(i)+1:length(logNROISum{i})),'.','Color','b'); %'+k' original symbol
plot(fullxsA,fg(coefsAy(:,i),fullxsA),'Color','r');
plot(fullxsB,fg(coefsBy(:,i),fullxsB),'Color','b');
hold off;
title({'979 Gauss',['Red a = ' num2str(coefsAy(1,i),3) '. Blue a = ' num2str(coefsBy(1,i),3) '.']});












