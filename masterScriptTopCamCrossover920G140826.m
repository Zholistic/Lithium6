directory = 'C:\Data\140822_crossover_920G_Isat0p3_10usPulse_freq5p4kHz_insitu\';
date = '140822';
camera = 'top';
varstring = 'motfet';
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
%for i=1:length(fileLocList)
%    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
%    imageArray(:,:,i) = atom2Image(:,:,1);
%end

calibrationImages = 30;

disp('Building images...');
for i=1:length(fileLocList(:))
    
    Isat = 10*135;
    
    %High Intensity images different isat:
    if(i > length(fileLocList(:))-calibrationImages) %12 Isat, high intensity images
        Isat = 135;  
    end
    
    %disp(num2str(Isat));
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArray(:,:,i) = atom2Image(:,:,1);
    
end

varDataLowIntensity = []; varDataHighIntensity = [];
varDataLowIntensity = varData(1:end-calibrationImages,1);
varDataHighIntensity = varData(end-calibrationImages+1:end,1);

%Crop images:
imageArrayC = []; imageArrayTC = []; imageArrayHighIntensityC = [];
imageArrayHighIntensityTC = [];
ROIx = 1:length(imageArray(1,:,1));
ROIy = 1:length(imageArray(:,1,1));
CrossROIy = 10:150; 
CrossROIx = 30:185; %The cross is inside the region specified above.
TightROIx = 30:185;
TightROIy = 10:150;
%Split into high and low intensity arrays
imageArrayC = imageArray(ROIy,ROIx,1:end-calibrationImages);
imageArrayHighIntensityC = imageArray(ROIy,ROIx,end-calibrationImages+1:end); 
imageArrayTC = imageArray(TightROIy,TightROIx,1:end-calibrationImages);
imageArrayHighIntensityTC = imageArray(TightROIy,TightROIx,end-calibrationImages+1:end); 

%Atom Number correction:
highIntImage = []; lowIntImage = [];
%There are 3 high intensity images for this dataset:
%0.77motfet:
highIntImage(:,:,1) = centerAndAverage(imageArrayHighIntensityC(:,:,1:10));
lowIntImage(:,:,1) = centerAndAverage(imageArrayC(:,:,427:431));
%1.01motfet:
highIntImage(:,:,2) = centerAndAverage(imageArrayHighIntensityC(:,:,11:20));
lowIntImage(:,:,2) = centerAndAverage(imageArrayC(:,:,367:371));
%1.23motfet:
highIntImage(:,:,3) = centerAndAverage(imageArrayHighIntensityC(:,:,21:end));
lowIntImage(:,:,3) = centerAndAverage(imageArrayC(:,:,206:220));

scaleFactor = [];
spectrumFunc1 = []; spectrumFunc2 = []; spectrumFunc3 = [];
spectrumFunc1 = makeSpectrumHL(highIntImage(:,:,1),lowIntImage(:,:,1));
spectrumFunc2 = makeSpectrumHL(highIntImage(:,:,2),lowIntImage(:,:,2));
spectrumFunc3 = makeSpectrumHL(highIntImage(:,:,3),lowIntImage(:,:,3));

scaleFactor(1) = mean(spectrumFunc1(1,50:82));
scaleFactor(2) = mean(spectrumFunc2(1,55:70));
scaleFactor(3) = mean(spectrumFunc3(1,60:70));

finalScaleFactor = mean(scaleFactor);
finalScaleFactor = 1.25; %No scaling due to on difference on these images (!?)

for i=1:length(imageArrayC(1,1,:))
    imageArrayC(:,:,i) = imageArrayC(:,:,i).*finalScaleFactor;
end


%Radially averaged profiles:
radProfiles = []; radProfilesT = []; center = [];
disp('Radially averaging...');
for i=1:length(imageArrayC(1,1,:))
    [radProfilesT(:,:,i),center(:,i)] = radAverageBigSquare(imageArrayC(:,:,i));
    radProfiles(:,:,i) = radProfilesT(:,1:end-20,i);
end

%Display every X image:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,5) == 0)       
        figure(i);
        imagesc(imageArrayC(:,:,i));        
    end
end
end


%%%%%Fit Gaussians:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
fgr = @(p,x)(p(1).*exp((-1).*((x).^2) ./ (2.*p(2).^2)));
gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
sigmaX = []; sigmaY = []; shiftFactor = []; shiftFactorR = []; gcoefsR = [];
disp('Gaussian Fitting...');
for i=1:length(imageArrayC(1,1,:))  
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    gcoefsR(:,i) = gausFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    
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
    sigmaR(:,i) = gcoefsR(2,i);
end

%Display every X fit:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsX(:,i),1:200)); hold on; plot(mean(imageArrayC(CrossROIy,:,i),1),'r'); hold off;      
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsY(:,i),1:180)); hold on; plot(mean(imageArrayC(:,CrossROIx,i),2),'r'); hold off;      
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fgr(gcoefsR(:,i),1:180)); hold on; plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); hold off;      
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
sortedVarData = []; indexs = []; sigmaRSort = [];
sigmaXSort = []; sigmaYSort = [];
[sortedVarData,indexs] = sort(varDataLowIntensity);
%indexs(:,1) is a vector of the sort.

for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
    sigmaRSort(i) = sigmaR(indexs(i));
end

%Average over same motfet data points:
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
widthsR = []; stdDevWidthsR = [];
sMomentXStdDev = []; sMomentYStdDev = [];
for i=1:length(sortedVarData)
    curr = sortedVarData(i);
    
    if( curr == prev )
        runTotal = runTotal+1;
    else
        %hit next value
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        widthsR(j) = mean(sigmaRSort(i-runTotal:i));
        stdDevWidthsR(j) = std(sigmaRSort(i-runTotal:i));
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
widthsX = widthsX.*2; %*2 to make it not the radius
widthsY = widthsY.*2;
widthsR = widthsR.*2;
%stdDevWidthsY = stdDevWidthsY.*pixelLength.*2; %full error on width
%stdDevWidthsX = stdDevWidthsX.*pixelLength.*2;
stdDevWidthsY = stdDevWidthsY.*2; %full error on width
stdDevWidthsX = stdDevWidthsX.*2;
stdDevWidthsR = stdDevWidthsR.*2;

%{
figure(1);
errorbar(pixelNumbers,widthsY,stdDevWidthsY/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('Y Widths vs Atom Number');
figure(2);
errorbar(pixelNumbers,widthsX,stdDevWidthsX/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('X Widths vs Atom Number');
%}
figure(3);
errorbar(motFets,pixelNumbers,pixelNumbersStdDev/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('Motfet vs Atom Number');
figure(4);
plot(log(pixelNumbers),log(widthsY),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('log(Y) Widths vs log(Atom Number)');
figure(5);
plot(log(pixelNumbers),log(widthsX),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('log(X) Widths vs log(Atom Number)');
figure(6);
errorbar(pixelNumbers,widthsY,stdDevWidthsY/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('Y Widths vs Atom Number');
hold on; plot(pixelCounts,gcoefsY(3,:)*2,'.r'); hold off;
figure(7);
errorbar(pixelNumbers,widthsX,stdDevWidthsX/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('X Widths vs Atom Number');
hold on; plot(pixelCounts,gcoefsX(3,:)*2,'.r'); hold off;
figure(8);
errorbar(pixelNumbers,widthsR,stdDevWidthsR/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('R Widths vs Atom Number');
hold on; plot(pixelCounts,gcoefsR(2,:).*2,'.r'); hold off;
figure(9);
plot(log(pixelNumbers),log(widthsR),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('Log R Widths vs Log Atom Number');
figure(10);
plot(pixelNumbers,widthsR./(pixelNumbers.^0.25),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('R Widths vs Atom Number');





