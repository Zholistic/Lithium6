directory = 'C:\Data\140911_2d_transverse_crossover_880G_750ms_laser_ramp_ISAT1e6_alpha_1\';
date = '140911';
camera = 'sidecam';
varstring = 'motfet';
magfield = '880G';
bins = 50;
%varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %2.84 um topcam, 
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
kB = 1.38e-23; %Boltzmanns constant m^2 kg s^-2 K^-1
PixelArea =(2.84e-6)^2;
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

varDataMain = [];
varDataMain = varData(:,1)';

imageArray = [];
%Pull images:
for i=1:length(fileLocList)
    imageArray(:,:,i) = PullFTS(fileLocList{i},raw);
end

%Crop images:
imageArrayC = []; imageArrayTC = [];
ROIx = 400:780;
ROIy = 500:675;
CrossROIx = 80:290; %The cross is inside the region specified above.
CrossROIy = 65:100;
TightROIx = 510:670;
TightROIy = 575:595;
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

%%%%%Fit Gaussians:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
fglz = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)));
%polylog order 1 function:
fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*x.^2)./(p(3).^2))));
gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
sigmaX = []; sigmaY = []; shiftFactor = []; halfXs = [];
disp('Function Fitting...');
for i=1:length(imageArrayC(1,1,:)) 
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    
    %Shift wings to zero:
    shiftFactor(i) = (gcoefsXi(4,i)+gcoefsYi(4,i))/2;
    imageArrayC(:,:,i) = imageArrayC(:,:,i) - shiftFactor(i);
    
    %Refit:
    gcoefsX(:,i) = gausFit1DLockZero(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    gcoefsY(:,i) = gausFit1DLockZero(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))
    
    centers(:,i) = [ceil(gcoefsX(2,i)), ceil(gcoefsY(2,i))]; %center = [x y]
    
    halfXs{i} = 1:(length(imageArrayC(1,:,i))-centers(1,i)+1);
    gcoefsPolyLog1(:,i) = polyLog1FitHalf1D(mean(imageArrayC(CrossROIy,centers(1,i):end,i),1),halfXs{i},camera);
    
    
    sigmaX(:,i) = gcoefsX(3,i);
    sigmaY(:,i) = gcoefsY(3,i);
end

%%%%%Second moment of fit function:
rvector = [];
rvector = halfXs;
COMr = 0; %Center of mass is located at zero
%for i=1:length(imageArrayTC(1,1,:))
%        %Sum of Intensity*pixel location / sum of intensity
%        COMx(i) = sum(mean(imageArrayTC(:,:,i),1).*xvector) / sum(mean(imageArrayTC(:,:,i),1));
%end

%Second moment direction
SMomR = []; sigmaRsm = [];
for i=1:length(imageArrayC(1,1,:))
    SMomR(i) = sum(fgp(gcoefsPolyLog1(:,i),rvector{i}).*(rvector{i} - COMr).^2) / sum(fgp(gcoefsPolyLog1(:,i),rvector{i}));
end

sigmaRsm = sqrt(SMomR);

%%%%%Temperatures:
TonTFs = [];
for i=1:length(imageArrayC(1,1,:))
    TonTFs(i) = 1/(log(1+exp(gcoefsPolyLog1(2,i)/gcoefsPolyLog1(3,i)^2)));
end

%Display every X fit:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsX(:,i),1:400));hold on; plot(mean(imageArrayC(CrossROIy,:,i),1),'r'); hold off;      
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

%Pixel to Atom Correction:
pixelCounts(:) = varData(:,5)*0.42;


%%%%%Culling of data points:
keepMe = []; j = 1;
for i=1:length(imageArrayC(1,1,:))
    if((pixelCounts(i) < 12000 && sigmaY(i) > 5) || (pixelCounts(i) < 12000 && sigmaY(i) < 4.1) )
        disp(['Discarded point number ' num2str(i)]);
    else
        keepMe(j) = i;
        j = j+1;
    end    
end

for i=1:length(keepMe)
    imageArrayC(:,:,i) = imageArrayC(:,:,keepMe(i));
    sigmaX(i) = sigmaX(keepMe(i));
    sigmaY(i) = sigmaY(keepMe(i));
    pixelCounts(i) = pixelCounts(keepMe(i));
    sigmaRsm(i) = sigmaRsm(keepMe(i));
    TonTFs(i) = TonTFs(keepMe(i));
end

%Sort varData:
indexs = [];
[sortedVarData,indexs] = sort(varDataMain);
%indexs(:,1) is a vector of the sort.
pixelCountsSort = []; sigmaXSort = []; sigmaYSort = [];
TonTFsSort = []; imageArrayCSort = [];
sigmaRsmSort = []; TonTFsSort = [];

for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
    
    sigmaRsmSort(i) = sigmaRsm(indexs(i));
    TonTFsSort(i) = TonTFs(indexs(i));
    imageArrayCSort(:,:,i) = imageArrayC(:,:,indexs(i));  
end

%Average over same motfet data points:
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
sMomentXStdDev = []; sMomentYStdDev = []; imageArrayAvgs = [];
widthsPsm = []; stdDevWidthsPsm = []; TonTFsm = []; stdDevTonTFsm = [];
prev = sortedVarData(1);
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
        motFets(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:i));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        runTotal = 0;
        j = j+1;
    end
    if(i == length(sortedVarData))
        %final iteration
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        motFets(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:i));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        runTotal = 0;
        j = j+1;
    end
    
    prev = curr;
end

if(0)
for i=1:length(imageArrayAvgs(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayAvgs(:,:,i));
    end
end
end


%Fit functions to the averaged images:
gcoefsXa = []; gcoefsYa = []; gcoefsYaError = []; gcoefsXaError = [];
for i=1:length(imageArrayAvgs(1,1,:))
    [gcoefsXa(:,i),gcoefsXaError(:,:,i)] = gausFit1D(mean(imageArrayAvgs(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    [gcoefsYa(:,i),gcoefsYaError(:,:,i)] = gausFit1D(mean(imageArrayAvgs(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))       
end

widthsYavg = gcoefsYa(3,:);
widthsXavg = gcoefsXa(3,:);

widthsYavgError = []; widthsXavgError = [];
for i=1:length(widthsYavg)
widthsYavgError(i) = widthsYavg(i) - gcoefsYaError(3,1,i);
widthsXavgError(i) = widthsXavg(i) - gcoefsXaError(3,1,i);
end

%Bin on atom number instead of motfet averaging:
%Raw arrays: sigmaX , sigmaY , pixelCounts, sigmaRsm, TonTFs, imageArrayC
sigmaXbinA = []; 
sigmaXbinA = binMeIncNaN(sigmaX,pixelCounts,bins);
% errorbar(sigmaXbinA(2,:),sigmaXbinA(1,:),sigmaXbinA(3,:),'.');
sigmaYbinA = [];
sigmaYbinA = binMeIncNaN(sigmaY,pixelCounts,bins);
% errorbar(sigmaYbinA(2,:),sigmaYbinA(1,:),sigmaYbinA(3,:),'.');
[avgImagesBin, atomNumsBin] = binMeCenterAndAverage(imageArrayC,pixelCounts,bins);

%Fit functions to the averaged bin images:
gcoefsXaBin = []; gcoefsYaBin = []; gcoefsYaBinError = []; gcoefsXaBinError = [];
for i=1:length(avgImagesBin(1,1,:))
    [gcoefsXaBin(:,i),gcoefsXaBinError(:,:,i)] = gausFit1DLockZero(mean(avgImagesBin(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    [gcoefsYaBin(:,i),gcoefsYaBinError(:,:,i)] = gausFit1DLockZero(mean(avgImagesBin(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))       
end

widthsYavgBin = gcoefsYaBin(3,:);
widthsXavgBin = gcoefsXaBin(3,:);

widthsYavgBinError = []; widthsXavgBinError = [];
for i=1:length(widthsYavgBin)
widthsYavgBinError(i) = widthsYavgBin(i) - gcoefsYaBinError(3,1,i);
widthsXavgBinError(i) = widthsXavgBin(i) - gcoefsXaBinError(3,1,i);
end

%convert to real units:
%widthsX = widthsX.*pixelLength.*2; %*2 to make it not the radius
%widthsY = widthsY.*pixelLength.*2; 
%stdDevWidthsY = stdDevWidthsY.*pixelLength.*2; %full error on width
%stdDevWidthsX = stdDevWidthsX.*pixelLength.*2;

figure(1);
errorbar(pixelNumbers,widthsY,stdDevWidthsY/2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCountsSort,sigmaYSort,'.r'); hold off;
figure(2);
errorbar(pixelNumbers,widthsX,stdDevWidthsX/2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCountsSort,sigmaXSort,'.r'); hold off;
figure(3);
errorbar(motFets,pixelNumbers,pixelNumbersStdDev/2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
figure(4);
plot(log(pixelNumbers),log(widthsX),'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
figure(13);
errorbar(pixelNumbers,widthsPsm,stdDevWidthsPsm./2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
title('Radial Widths vs Pixel Number, Second Moment of PolyLog Fit');
figure(16);
errorbar(pixelNumbers,TonTFsm,stdDevTonTFsm./2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
title('TonTF vs Pixel Number');
figure(20)
errorbar(pixelNumbers,widthsYavg,widthsYavgError,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCountsSort,sigmaYSort,'.r'); hold off;
figure(21)
errorbar(pixelNumbers,widthsXavg,widthsXavgError,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
figure(25)
errorbar(atomNumsBin,widthsYavgBin,widthsYavgBinError,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCountsSort,sigmaYSort,'.r'); hold off;
figure(26)
errorbar(atomNumsBin,widthsXavgBin,widthsXavgBinError,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
figure(30);
errorbar(sigmaXbinA(2,:),sigmaXbinA(1,:),sigmaXbinA(3,:),'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCounts,sigmaX,'.r'); hold off;
h = figure(31);
errorbar(sigmaYbinA(2,:),sigmaYbinA(1,:),sigmaYbinA(3,:),'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCountsSort,sigmaYSort,'.r'); hold off;
figname = [date '_' camera '_' magfield '_Tight_' num2str(bins) 'Bins'];
figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\';
saveas(h,[figdirectory figname '.fig'],'fig');
saveas(h,[figdirectory figname '.png'],'png');

if(0)
%Binning Debug!:
figure(50);
errorbar(sigmaYbinA(2,:),sigmaYbinA(1,:),sigmaYbinA(3,:),'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
hold on; plot(pixelCountsSort,sigmaYSort,'.r'); 
for i=1:length(sigmaYbinA(4,:))
line([sigmaYbinA(4,i) sigmaYbinA(4,i)],[min(sigmaYbinA(1,:)) max(sigmaYbinA(1,:))],'LineStyle','--','Color',[0.7 0.7 0.7]); 
end
%for i=1:length(binEdge)
%line([binEdge(i) binEdge(i)],[min(sigmaYbinA(1,:)) max(sigmaYbinA(1,:))],'LineStyle','--','Color',[0.1 0.8 0.1]); 
%end
hold off;
end


