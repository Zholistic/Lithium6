directory = 'C:\Data\140818_crossover_972G_Isat0p2_9p44usPulse_freq5p4kHz_insitu\';
date = '140818';
camera = 'top';
varstring = 'motfet ';
%varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %2.84 um topcam, 
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
kB = 1.38e-23; %Boltzmanns constant m^2 kg s^-2 K^-1
imgArrayFresh = [];  lowIntRealAtomImg = [];
OD = 0; %optical density from SPE process function 1=OD, 0=WithSigma
close all;
raw = 1;


%Read in the log file:

logfilename = [directory date '_log_camera.txt'];
fid = fopen(logfilename,'rt');
C = textscan(fid, '%s', 'Delimiter','\t'); %tokenize into tab seperated tokens
C = C{1};
fclose(fid);

%Iterate through strings in the log file array 1-d C{:} to find the
%information we want:

varData = [];

TotalImages = 0;

%for loop to find the total number of images
for i=1:(length(C)-1)

    curr = C{i};
    next = C{i+1};
    
    if strcmp(curr,'MeasNr')
        TotalImages = TotalImages + 1;
    end
              
    
end

%Get information from log file:
[fileLocList,varData] = generateFromLogfile(directory,date,varstring,camera);

imageArray = [];
%Pull images:
for i=1:length(fileLocList)
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArray(:,:,i) = atom2Image(:,:,1);
end

%Crop images:
imageArrayC = []; imageArrayTC = [];
ROIx = 1:length(imageArray(1,:,1));
ROIy = 1:length(imageArray(:,1,1));
CrossROIy = 10:150; 
CrossROIx = 30:185; %The cross is inside the region specified above.
TightROIx = 30:185;
TightROIy = 10:150;
imageArrayC = imageArray(ROIy,ROIx,:);
imageArrayTC = imageArray(TightROIy,TightROIx,:);


%Display every X image:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,5) == 0)       
        figure(i);
        imagesc(imageArrayC(:,:,i));        
    end
end
end

%Radially Average:
%for i=1:length(imageArrayC)
%radProfile = radAverageBigSquare(lowIntRealAtomImg);
%figure(203)
%plot(radProfile(2,:),radProfile(1,:));

%%%%%Fit Gaussians:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
sigmaX = []; sigmaY = []; shiftFactor = [];
disp('Gaussian Fitting...');
for i=1:length(imageArrayC(1,1,:))  
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    
    %Shift wings to zero:
    shiftFactor(i) = (gcoefsXi(4,i)+gcoefsYi(4,i))/2;
    imageArrayC(:,:,i) = imageArrayC(:,:,i) - shiftFactor(i);
    
    %Refit:
    gcoefsX(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    gcoefsY(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))
    
    centers(:,i) = [ceil(gcoefsX(2,i)), ceil(gcoefsY(2,i))]; %center = [x y]
    sigmaX(:,i) = gcoefsX(3,i);
    sigmaY(:,i) = gcoefsY(3,i);
end

%Display every X fit:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsX(:,i),1:200));hold on; plot(mean(imageArrayC(CrossROIy,:,i),1),'r'); hold off;      
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsY(:,i),1:180));hold on; plot(mean(imageArrayC(:,CrossROIx,i),2),'r'); hold off;      
    end
end
end

%%%%%Atom numbers:
%Tight ROI array:
pixelCounts = [];
for i=1:length(imageArrayC(1,1,:))
    pixelCounts(i) = sum(sum(imageArrayC(:,:,i)));
end


%Sort varData:
[sortedVarData,indexs] = sort(varData);
%indexs(:,1) is a vector of the sort.

for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
end

%Average over same motfet data points:
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
sMomentXStdDev = []; sMomentYStdDev = [];
for i=1:length(sortedVarData)
    curr = sortedVarData(i,1);
    
    if( curr == prev )
        runTotal = runTotal+1;
    else
        %hit next value
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        motFets(j) = sortedVarData(i);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        runTotal = 0;
        j = j+1;
    end
    
    prev = curr;
end

%convert to real units:
%widthsX = widthsX.*pixelLength.*2; %*2 to make it not the radius
%widthsY = widthsY.*pixelLength.*2; 
%stdDevWidthsY = stdDevWidthsY.*pixelLength.*2; %full error on width
%stdDevWidthsX = stdDevWidthsX.*pixelLength.*2;

figure(1);
errorbar(pixelNumbers,widthsY,stdDevWidthsY/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(2);
errorbar(pixelNumbers,widthsX,stdDevWidthsX/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(3);
errorbar(motFets,pixelNumbers,pixelNumbersStdDev/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(4);
plot(log(pixelNumbers),log(widthsX),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;




