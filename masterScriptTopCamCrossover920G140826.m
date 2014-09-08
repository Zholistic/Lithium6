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
lowIntImage(:,:,1) = centerAndAverage(imageArrayC(:,:,426:430));
%1.01motfet:
highIntImage(:,:,2) = centerAndAverage(imageArrayHighIntensityC(:,:,11:20));
lowIntImage(:,:,2) = centerAndAverage(imageArrayC(:,:,366:370));
%1.23motfet:
highIntImage(:,:,3) = centerAndAverage(imageArrayHighIntensityC(:,:,21:end));
lowIntImage(:,:,3) = centerAndAverage(imageArrayC(:,:,205:219));

scaleFactor = [];
spectrumFunc1 = []; spectrumFunc2 = []; spectrumFunc3 = [];
spectrumFunc1 = makeSpectrumHL(highIntImage(:,:,1),lowIntImage(:,:,1));
spectrumFunc2 = makeSpectrumHL(highIntImage(:,:,2),lowIntImage(:,:,2));
spectrumFunc3 = makeSpectrumHL(highIntImage(:,:,3),lowIntImage(:,:,3));

scaleFactor(1) = mean(spectrumFunc1(1,50:82));
scaleFactor(2) = mean(spectrumFunc2(1,55:70));
scaleFactor(3) = mean(spectrumFunc3(1,60:70));

finalScaleFactor = mean(scaleFactor);
finalScaleFactor = 1.35; %No scaling due to on difference on these images (!?)

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

%Shift wings of radial profiles to zero:
shiftBy = [];
for i=1:length(radProfiles(1,1,:))
    shiftBy(i) = mean(radProfiles(1,55:75,i));
    radProfiles(1,:,i) = radProfiles(1,:,i) - shiftBy(i);
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


%%%%%Fit Functions:
%gaussian with 4 variables:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
%gaussian with 2 variables (for radial profiles)
fgr = @(p,x)(p(1).*exp((-1).*((x).^2) ./ (2.*p(2).^2)));
%polylog order 1 function:
fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*x.^2)./(p(3).^2))));
%gaussian with 2 variables fixed 2.5 exponent:
fgr2p5 = @(p,x)(p(1).*exp((-1).*((x).^(2.5)) ./ (2.*p(2).^(2.5))));
%gaussian with 2 variables and fitted exponent:
fgrfe = @(p,x)(p(1).*exp((-1).*((x).^(p(2))) ./ (2.*p(3).^(2.5)))); 

gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
gcoefsR = []; gcoefsR2p5 = []; sigmaR2p5 = []; gcoefsPolyLog1 = [];
gcoefsRFreeExp = [];
sigmaX = []; sigmaY = []; shiftFactor = []; shiftFactorR = [];
disp('Function Fitting...');
for i=1:length(imageArrayC(1,1,:))
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    gcoefsR(:,i) = gausFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    %gcoefsR2p5(:,i) = gausExp2p5FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    gcoefsPolyLog1(:,i) = polyLog1FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    %gcoefsRFreeExp(:,i) = gausFreeExpFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    
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
    %sigmaR2p5(:,i) = gcoefsR2p5(2,i);
end


%%%%%Second moment of fit function:
rvector = [];
rvector = radProfiles(2,:,1);
COMr = 0; %Center of mass is located at zero
%for i=1:length(imageArrayTC(1,1,:))
%        %Sum of Intensity*pixel location / sum of intensity
%        COMx(i) = sum(mean(imageArrayTC(:,:,i),1).*xvector) / sum(mean(imageArrayTC(:,:,i),1));
%end

%Second moment direction
SMomR = []; sigmaRsm = [];
for i=1:length(radProfiles(1,1,:))
    SMomR(i) = sum(fgp(gcoefsPolyLog1(:,i),rvector).*(rvector - COMr).^2) / sum(fgp(gcoefsPolyLog1(:,i),rvector));
end

sigmaRsm = sqrt(SMomR);


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
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,6) == 0)       
        figure(i);
        plot(fgr2p5(gcoefsR2p5(:,i),1:180),'g'); hold on; plot(fgr(gcoefsR(:,i),1:180)); plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); hold off;      
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,6) == 0)       
        figure(i);
        plot(fgr(gcoefsR(:,i),1:180),'g'); hold on; plot(fgp(gcoefsPolyLog1(:,i),1:180)); plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); line([sigmaRsm(i) sigmaRsm(i)],[0 fgp(gcoefsPolyLog1(:,i),sigmaRsm(i))],'LineStyle','--','Color',[0.7 0.7 0.7]); hold off;      
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,6) == 0)       
        figure(i);
        plot(fgr(gcoefsR(:,i),1:180),'g'); hold on; plot(fgrfe(gcoefsRFreeExp(:,i),1:180)); plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); hold off;      
    end
end
end

%%%%%Atom numbers:
%Tight ROI array:
pixelCounts = [];
for i=1:length(imageArrayC(1,1,:))
    pixelCounts(i) = sum(sum(imageArrayC(:,:,i)));
end

%%%%%Temperatures:
TonTFs = [];
for i=1:length(imageArrayC(1,1,:))
    TonTFs(i) = 1/(log(1+exp(gcoefsPolyLog1(2,i)/gcoefsPolyLog1(3,i)^2)));
end

%Sort varData:
sortedVarData = []; indexs = []; sigmaRSort = [];
sigmaXSort = []; sigmaYSort = []; sigmaR2p5Sort = [];
sigmaRsmSort = []; TonTFsSort = []; imageArrayCSort = [];
[sortedVarData,indexs] = sort(varDataLowIntensity);
%indexs(:,1) is a vector of the sort.

for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
    sigmaRSort(i) = sigmaR(indexs(i));
    %sigmaR2p5Sort(i) = sigmaR2p5(indexs(i));
    sigmaRsmSort(i) = sigmaRsm(indexs(i));
    TonTFsSort(i) = TonTFs(indexs(i));
    imageArrayCSort(:,:,i) = imageArrayC(:,:,indexs(i));  
end

%Average over same motfet data points:
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
widthsR = []; stdDevWidthsR = []; widthsR2p5 = [];
widthsPsm = []; stdDevWidthsPsm = []; TonTFsm = []; stdDevTonTFsm = [];
sMomentXStdDev = []; sMomentYStdDev = []; stdDevWidthsR2p5 = [];
prev = sortedVarData(1); imageArrayAvgs = [];
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
        %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
        %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
        widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:i));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        motFets(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        runTotal = 0;
        j = j+1;
    end
    if( i == length(sortedVarData))
        %Final run:
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        widthsR(j) = mean(sigmaRSort(i-runTotal:i));
        stdDevWidthsR(j) = std(sigmaRSort(i-runTotal:i));
        %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
        %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
        widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:i));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        motFets(j) = sortedVarData(i-runTotal);
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
%widthsR2p5 = widthsR2p5.*2;
widthsPsm = widthsPsm.*2;
%stdDevWidthsY = stdDevWidthsY.*pixelLength.*2; %full error on width
%stdDevWidthsX = stdDevWidthsX.*pixelLength.*2;
stdDevWidthsY = stdDevWidthsY.*2; %full error on width
stdDevWidthsX = stdDevWidthsX.*2;
stdDevWidthsR = stdDevWidthsR.*2;
%stdDevWidthsR2p5 = stdDevWidthsR2p5.*2;
stdDevWidthsPsm = stdDevWidthsPsm.*2;

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
if(0)
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
errorbar(pixelNumbers,widthsY,stdDevWidthsY./2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('Y Widths vs Atom Number');
hold on; plot(pixelCounts,gcoefsY(3,:)*2,'.r'); hold off;
figure(7);
errorbar(pixelNumbers,widthsX,stdDevWidthsX./2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('X Widths vs Atom Number');
hold on; plot(pixelCounts,gcoefsX(3,:)*2,'.r'); hold off;
figure(8);
errorbar(pixelNumbers,widthsR,stdDevWidthsR./2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
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
figure(11);
plot(pixelNumbers,widthsR./(pixelNumbers.^0.25),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('R Widths/N^{0.25} vs Atom Number');
end
%figure(12);
%errorbar(pixelNumbers,widthsR2p5,stdDevWidthsR2p5./2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
%    'Marker','o',...
%    'LineStyle','none',...
%    'Color',[0 0 1]);
%grid on;
%title('R Widths vs Atom Number, 2.5 exponent');
figure(13);
errorbar(pixelNumbers,widthsPsm,stdDevWidthsPsm./2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('R Widths vs Atom Number, Second Moment of PolyLog Fit');
hold on; plot(pixelCounts,sigmaRsm.*2,'.r'); hold off;
figure(14);
plot(log10(pixelNumbers),log10(widthsPsm),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('Log10 R Widths vs Log10 Atom Number, Second Moment of PolyLog Fit');
figure(15);
plot(pixelNumbers,widthsPsm./(pixelNumbers.^0.25),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCounts,(sigmaRsm.*2)./(pixelCounts.^(1/4)),'.r'); hold off;
title('R Widths/N^{0.25} vs Atom Number, Second Moment of PolyLog Fit');
figure(16);
errorbar(pixelNumbers,TonTFsm,stdDevTonTFsm./2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
title('TonTF vs Atom Number');


