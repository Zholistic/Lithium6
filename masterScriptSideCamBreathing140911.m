
directory = 'C:\Data\140911_breathingvsatomnumber_50msRamp\';
date = '140911';
camera = 'top';
varstring = 'holdtime';
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

for i=1:length(fileLocList(:))
    
    Isat = 10*135;
    
    %High Intensity images different isat:
    if(0) %NO High intensity data for this run
        Isat = 135;  
    end
    
    %disp(num2str(Isat));
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArray(:,:,i) = atom2Image(:,:,1);
    
end

varDataLowIntensity = []; varDataHighIntensity = [];
varDataLowIntensity = varData(:,1)';
%varDataHighIntensity = varData(381:end)';

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
imageArrayC = imageArray(ROIy,ROIx,:);
%imageArrayHighIntensityC = imageArray(ROIy,ROIx,381:end); 
imageArrayTC = imageArray(TightROIy,TightROIx,:);
%imageArrayHighIntensityTC = imageArray(TightROIy,TightROIx,381:end); 

if(0)
%Atom Number correction:
highIntImage = []; lowIntImage = [];
%There are 0 high intensity images for this dataset:
%0.89motfet:
highIntImage(:,:,1) = centerAndAverage(imageArrayHighIntensityC(:,:,1:10));
lowIntImage(:,:,1) = centerAndAverage(imageArrayC(:,:,332:336));
%1.1motfet:
highIntImage(:,:,2) = centerAndAverage(imageArrayHighIntensityC(:,:,10:end));
lowIntImage(:,:,2) = centerAndAverage(imageArrayC(:,:,257:277));

scaleFactor = [];
spectrumFunc1 = []; spectrumFunc2 = [];
spectrumFunc1 = makeSpectrumHL(highIntImage(:,:,1),lowIntImage(:,:,1));
spectrumFunc2 = makeSpectrumHL(highIntImage(:,:,2),lowIntImage(:,:,2));

scaleFactor(1) = mean(spectrumFunc1(1,45:75));
scaleFactor(2) = mean(spectrumFunc2(1,40:85));

finalScaleFactor = mean(scaleFactor);
end

finalScaleFactor = 1.3;

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
    gcoefsPolyLog1(:,i) = polyLog1FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i),camera);
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

%Datasets:
dataset1 = 1:44;
dataset2 = 45:82;
dataset3 = 84:104;
dataset4 = 105:125;

holdtimes = [];
holdtimes{1} = varData(dataset1,1);
holdtimes{2} = varData(dataset2,1);
holdtimes{3} = varData(dataset3,1);
holdtimes{4} = varData(dataset4,1);

sigmaRsmSort = [];
sigmaRsmSort = sigmaRsm;

if(0)
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

end

%Average over same data points:
holdTimes = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
widthsR = []; stdDevWidthsR = []; widthsR2p5 = [];
widthsPsm = []; stdDevWidthsPsm = []; TonTFsm = []; stdDevTonTFsm = [];
sMomentXStdDev = []; sMomentYStdDev = []; stdDevWidthsR2p5 = [];
imageArrayAvgs = [];

for set=1:4
    
    sortedVarData = holdtimes{set};    
    j=1; runTotal = 0; prev = sortedVarData(1);
    
    for i=1:length(sortedVarData)
        curr = sortedVarData(i);
        
        if( curr == prev )
            runTotal = runTotal+1;
        elseif (j == 1)
            %special case for if there is only one datapoint
            %in the first set:
            widthsX(j,set) = mean(sigmaXSort(i-runTotal:j));
            stdDevWidthsX(j,set) = std(sigmaXSort(i-runTotal:j));
            widthsY(j,set) = mean(sigmaYSort(i-runTotal:j));
            stdDevWidthsY(j,set) = std(sigmaYSort(i-runTotal:j));
            widthsR(j,set) = mean(sigmaRSort(i-runTotal:j));
            stdDevWidthsR(j,set) = std(sigmaRSort(i-runTotal:j));
            %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
            %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
            widthsPsm(j,set) = mean(sigmaRsmSort(i-runTotal:j));
            stdDevWidthsPsm(j,set) = std(sigmaRsmSort(i-runTotal:j));
            
            TonTFsm(j,set) = mean(TonTFsSort(i-runTotal:j));
            stdDevTonTFsm(j,set) = std(TonTFsSort(i-runTotal:j));
            
            %imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
            
            holdTimes(j,set) = sortedVarData(i-runTotal);
            pixelNumbers(j,set) = mean(pixelCountsSort(i-runTotal:j));
            pixelNumbersStdDev(j,set) = std(pixelCountsSort(i-runTotal:j));
            runTotal = 0;
            j = j+1;
        else
            %hit next value
            widthsX(j,set) = mean(sigmaXSort(i-runTotal:i));
            stdDevWidthsX(j,set) = std(sigmaXSort(i-runTotal:i));
            widthsY(j,set) = mean(sigmaYSort(i-runTotal:i));
            stdDevWidthsY(j,set) = std(sigmaYSort(i-runTotal:i));
            widthsR(j,set) = mean(sigmaRSort(i-runTotal:i));
            stdDevWidthsR(j,set) = std(sigmaRSort(i-runTotal:i));
            %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
            %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
            widthsPsm(j,set) = mean(sigmaRsmSort(i-runTotal:i));
            stdDevWidthsPsm(j,set) = std(sigmaRsmSort(i-runTotal:i));
            
            TonTFsm(j,set) = mean(TonTFsSort(i-runTotal:i));
            stdDevTonTFsm(j,set) = std(TonTFsSort(i-runTotal:i));
            
            %imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
            
            holdTimes(j,set) = sortedVarData(i-runTotal);
            pixelNumbers(j,set) = mean(pixelCountsSort(i-runTotal:i));
            pixelNumbersStdDev(j,set) = std(pixelCountsSort(i-runTotal:i));
            runTotal = 0;
            j = j+1;
        end
        if( i == length(sortedVarData))
            %Final run:
            widthsX(j,set) = mean(sigmaXSort(i-runTotal:i));
            stdDevWidthsX(j,set) = std(sigmaXSort(i-runTotal:i));
            widthsY(j,set) = mean(sigmaYSort(i-runTotal:i));
            stdDevWidthsY(j,set) = std(sigmaYSort(i-runTotal:i));
            widthsR(j,set) = mean(sigmaRSort(i-runTotal:i));
            stdDevWidthsR(j,set) = std(sigmaRSort(i-runTotal:i));
            %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
            %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
            widthsPsm(j,set) = mean(sigmaRsmSort(i-runTotal:i));
            stdDevWidthsPsm(j,set) = std(sigmaRsmSort(i-runTotal:i));
            
            TonTFsm(j,set) = mean(TonTFsSort(i-runTotal:i));
            stdDevTonTFsm(j,set) = std(TonTFsSort(i-runTotal:i));
            
            %imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
            
            holdTimes(j,set) = sortedVarData(i-runTotal);
            pixelNumbers(j,set) = mean(pixelCountsSort(i-runTotal:i));
            pixelNumbersStdDev(j,set) = std(pixelCountsSort(i-runTotal:i));
            runTotal = 0;
            j = j+1;
        end
        
        prev = curr;
    end
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


figure(1);
errorbar(holdTimes(:,1),widthsPsm(:,1)./2,stdDevWidthsPsm(:,1)./4,'MarkerSize',5,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(2);
plot(varData(dataset2,1),sigmaRsm(dataset2),'MarkerSize',5,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(3);
plot(varData(dataset3,1),sigmaRsm(dataset3),'MarkerSize',5,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(4);
plot(varData(dataset4,1),sigmaRsm(dataset4),'MarkerSize',5,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;






