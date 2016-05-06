directory = 'C:\Data\160505_CollectiveOscillations_2D_Green_modulation_triangle_envelope_880G_20_oscillation\';

date = '160505';
camera = 'topcam';
varstring = 'HoldTime';
%varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 134*10; %135*x us
%Isat = 10^6;
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
widthYlogfile = []; widthXlogfile = []; nROIlogfile = [];
[fileLocList,varData] = generateFromLogfile(directory,date,varstring,camera);
widthYlogfile = varData(:,4);
widthXlogfile = varData(:,3);
nROIlogfile = varData(:,5);

imageArray = [];
%Pull images:
for i=1:length(fileLocList)
    OD = 0;
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArray(:,:,i) = atom2Image(:,:,1);
end

%Crop images:
imageArrayC = []; imageArrayTC = [];
ROIx = 1:length(imageArray(1,:,1));
ROIy = 1:length(imageArray(:,1,1));
CrossROIy = 20:180; 
CrossROIx = 35:185; %The cross is inside the region specified above.
TightROIx = 35:185;
TightROIy = 30:170;
%Split into high and low intensity arrays
imageArrayC = imageArray(ROIy,ROIx,:);
imageArrayTC = imageArray(TightROIy,TightROIx,:);

%Display every X image:
if(0)
for i=1:length(imageArray(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArray(:,:,i));        
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayTC(:,:,i));        
    end
end
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

%%%%%Fit Functions:
%gaussian with 4 variables:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
%gaussian with 2 variables (for radial profiles)
fgr = @(p,x)(p(1).*exp((-1).*((x).^2) ./ (2.*p(2).^2)));
%polylog order 1 function:
%fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*x.^2)./(p(3)^2))));
fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*p(3).*x.^2)./(p(4)))));
%gaussian with 2 variables fixed 2.5 exponent:
fgr2p5 = @(p,x)(p(1).*exp((-1).*((x).^(2.5)) ./ (2.*p(2).^(2.5))));
%gaussian with 2 variables and fitted exponent:
fgrfe = @(p,x)(p(1).*exp((-1).*((x).^(p(2))) ./ (2.*p(3).^(2.5)))); 

gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
gcoefsR = []; gcoefsR2p5 = []; sigmaR2p5 = []; gcoefsPolyLog1 = [];
gcoefsRFreeExp = []; sigmaR = [];
sigmaX = []; sigmaY = []; shiftFactor = []; shiftFactorR = [];
disp('Function Fitting...');
for i=1:length(imageArrayC(1,1,:))
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    gcoefsR(:,i) = gausFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    %gcoefsR2p5(:,i) = gausExp2p5FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    gcoefsPolyLog1(:,i) = polyLog1FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i),'topcam');
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


%show function fits:
if(0)
    for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        plot(fgr(gcoefsR(:,i),1:180),'g'); hold on; plot(fgp(gcoefsPolyLog1(:,i),1:180)); plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); line([sigmaRsm(i) sigmaRsm(i)],[0 fgp(gcoefsPolyLog1(:,i),sigmaRsm(i))],'LineStyle','--','Color',[0.7 0.7 0.7]); hold off;      
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
TonTFs = []; Ts = [];
for i=1:length(imageArrayC(1,1,:))
    TonTFs(i) = 1/(log(1+exp(gcoefsPolyLog1(2,i)/gcoefsPolyLog1(3,i)^2)));
    Ts(i) = gcoefsPolyLog1(4,i);
end


%Sort varData:
sortedVarData = []; indexs = []; sigmaRSort = [];
sigmaXSort = []; sigmaYSort = []; sigmaR2p5Sort = [];
sigmaRsmSort = []; TonTFsSort = []; imageArrayCSort = [];
[sortedVarData,indexs] = sort(varData);
%indexs(:,1) is a vector of the sort (varstring).

%Re-order each array:
for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
    sigmaRSort(i) = sigmaR(indexs(i));
    %sigmaR2p5Sort(i) = sigmaR2p5(indexs(i));
    %sigmaRsmSort(i) = sigmaRsm(indexs(i));
    TonTFsSort(i) = TonTFs(indexs(i));
    imageArrayCSort(:,:,i) = imageArrayC(:,:,indexs(i));  
end

%Average over same motfet data points:
%Bug fix history: 150109 fixed endNum
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
widthsR = []; stdDevWidthsR = []; widthsR2p5 = [];
widthsPsm = []; stdDevWidthsPsm = []; TonTFsm = []; stdDevTonTFsm = [];
sMomentXStdDev = []; sMomentYStdDev = []; stdDevWidthsR2p5 = [];
prev = sortedVarData(1); imageArrayAvgs = []; runTotals = [];
for i=1:length(sortedVarData)
    curr = sortedVarData(i);
    
    if( curr == prev )
        runTotal = runTotal+1;
    else
        %hit next value
        runTotals(j) = runTotal;
        startNum = i-runTotal;
        endNum = i-1;
        widthsX(j) = mean(sigmaXSort(i-runTotal:endNum));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:endNum));
        widthsY(j) = mean(sigmaYSort(i-runTotal:endNum));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:endNum));
        widthsR(j) = mean(sigmaRSort(i-runTotal:endNum));
        stdDevWidthsR(j) = std(sigmaRSort(i-runTotal:endNum));
        %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
        %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
        %widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        %stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:endNum));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:endNum));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:endNum));
        
        motFets(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:endNum));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:endNum));
        runTotal = 1;
        j = j+1;
    end
    if( i == length(sortedVarData))
        %Final run:
        if(i == runTotal)
            runTotal = runTotal-1;
        end
        runTotals(j) = runTotal;
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        widthsR(j) = mean(sigmaRSort(i-runTotal:i));
        stdDevWidthsR(j) = std(sigmaRSort(i-runTotal:i));
        %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
        %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
        %widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        %stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:i));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        motFets(j) = sortedVarData(i);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        runTotal = 1;
        j = j+1;
    end
    
    prev = curr;
end

if(0)
errorbar(pixelNumbers,widthsR,stdDevWidthsR./2,stdDevWidthsR./2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
end


if(0)
    for i=1:length(imageArrayAvgs(1,1,:))
    if(mod(i,2) == 0)       
        figure(i);
        imagesc(imageArrayAvgs(:,:,i));        
    end
    end
end


if(0)
    radProfileAvg = []; gcoefsRAvg = [];
    for i=1:length(imageArrayAvgs(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        
        [radProfilesAvg(:,:,i),~] = radAverageBigSquare(imageArrayAvgs(:,:,i));
        gcoefsRAvg(:,i) = gausFitHalf1D(radProfilesAvg(1,:,i),radProfilesAvg(2,:,i));
        plot(fgr(gcoefsRAvg(:,i),1:200));hold on; grid on; plot(radProfilesAvg(2,:,i), radProfilesAvg(1,:,i),'r'); hold off;      
    end
    end
end

plot(motFets(1:end),widthsR(1:end),'.');