directory = 'C:\Data\151005_breathing_2D_ramp_250ms_fieldramp_150ms_hold_150ms_compress_10ms_hold_5ms_decompress_0p1us\';
%directory = 'C:\Data\151003_New_breath_test\';
%directory = 'C:\Data\150904_magnification_check3_top\';
%NOT ACTUALLY 724G. like 780 or 750, check lab book...
date = '151005';
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
%Image 10 has the streaks
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

errorbar(pixelNumbers,widthsR,stdDevWidthsR./2,stdDevWidthsR./2,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;


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

%datasets: 1:70, 71:140, 141:210, 211:end
fgsine = @(p,x)(p(1).*sin(p(2).*x+p(3))+p(4));
fgsineDamp = @(p,x)(p(1).*exp(-p(2).*x).*sin(p(3).*x+p(4))+p(5));
if(0)
    set1 = 1:70; set2 = 71:140; set3 = 141:210; set4 = 211:280;
    datasx1 = varData(1:70,1); datasy1 = sigmaR(1:70); imagesSet1 = imageArrayC(:,:,1:70);
    datasx2 = varData(71:140,1); datasy2 = sigmaR(71:140); imagesSet2 = imageArrayC(:,:,71:140);
    datasx3 = varData(141:210,1); datasy3 = sigmaR(141:210); imagesSet3 = imageArrayC(:,:,141:210);
    datasx4 = varData(211:end,1); datasy4 = sigmaR(211:end); imagesSet4 = imageArrayC(:,:,211:end);
    [datasx1Sorted,indexsx1] = sort(datasx1);
    [datasx2Sorted,indexsx2] = sort(datasx2);
    [datasx3Sorted,indexsx3] = sort(datasx3);
    [datasx4Sorted,indexsx4] = sort(datasx4);
    for i=1:length(datasx1Sorted)
        datasy1Sorted(i) = datasy1(indexsx1(i));
        imagesSortedSet1(:,:,i) = imagesSet1(:,:,indexsx1(i));
    end
    for i=1:length(datasx2Sorted)
        datasy2Sorted(i) = datasy2(indexsx2(i));
        imagesSortedSet2(:,:,i) = imagesSet2(:,:,indexsx2(i));
    end
    for i=1:length(datasx3Sorted)
        datasy3Sorted(i) = datasy3(indexsx3(i));
        imagesSortedSet3(:,:,i) = imagesSet3(:,:,indexsx3(i));
    end
    for i=1:length(datasx4Sorted)
        datasy4Sorted(i) = datasy4(indexsx4(i));
        imagesSortedSet4(:,:,i) = imagesSet4(:,:,indexsx4(i));
    end
    
    %avgImagesBin = []; atomNumsBin = []; xElementsCell = [];
    %bins = ceil(length(datasx1)/2);
    %[avgImagesBin, atomNumsBin, xElementsCell] = binMeCenterAndAverageIncNaN(imageArrayC(:,:,1:70),pixelCounts(1:70),bins);
    %[radProfilesAvg(:,:,i),~] = radAverageBigSquare(imageArrayAvgs(:,:,i));
    %gcoefsRAvg(:,i) = gausFitHalf1D(radProfilesAvg(1,:,i),radProfilesAvg(2,:,i));
    
    %2 data points per timestep:
    j=1; datay1ssE = []; datay1ssEstdDev = []; datax1ssE = []; imagesAvgSet1 = [];
    for i=1:2:length(datasx1Sorted)
        datay1ssE(j) = mean(datasy1Sorted(i:i+1));
        datay1ssEstdDev(j) = std(datasy1Sorted(i:i+1));
        imagesAvgSet1(:,:,j) = centerAndAverage(imagesSortedSet1(:,:,i:i+1));
        datax1ssE(j) = datasx1Sorted(i);
        j = j+1;
        %i
    end
    
    j=1; datay2ssE = []; datay2ssEstdDev = []; datax2ssE = []; imagesAvgSet2 = [];
    for i=1:2:length(datasx2Sorted)
        datay2ssE(j) = mean(datasy2Sorted(i:i+1));
        datay2ssEstdDev(j) = std(datasy2Sorted(i:i+1));
        imagesAvgSet2(:,:,j) = centerAndAverage(imagesSortedSet2(:,:,i:i+1));
        datax2ssE(j) = datasx2Sorted(i);
        j = j+1;
        %i
    end
    
    j=1; datay3ssE = []; datay3ssEstdDev = []; datax3ssE = []; imagesAvgSet3 = [];
    for i=1:2:length(datasx3Sorted)
        datay3ssE(j) = mean(datasy3Sorted(i:i+1));
        datay3ssEstdDev(j) = std(datasy3Sorted(i:i+1));
        imagesAvgSet3(:,:,j) = centerAndAverage(imagesSortedSet3(:,:,i:i+1));
        datax3ssE(j) = datasx3Sorted(i);
        j = j+1;
        %i
    end
    
    j=1; datay4ssE = []; datay4ssEstdDev = []; datax4ssE = []; imagesAvgSet4 = [];
    for i=1:2:length(datasx4Sorted)
        datay4ssE(j) = mean(datasy4Sorted(i:i+1));
        datay4ssEstdDev(j) = std(datasy4Sorted(i:i+1));
        imagesAvgSet4(:,:,j) = centerAndAverage(imagesSortedSet4(:,:,i:i+1));
        datax4ssE(j) = datasx4Sorted(i);
        j = j+1;
        %i
    end
    
    %plot(datax1ssE,datay1ssE,'.'); hold on; plot(datax1ssE,datay1ssEstdDev);
    errorbar(datax1ssE,datay1ssE,datay1ssEstdDev./2,'.');
    errorbar(datax2ssE,datay2ssE,datay2ssEstdDev./2,'.r');
    errorbar(datax3ssE,datay3ssE,datay3ssEstdDev./2,'.g');
    errorbar(datax4ssE,datay4ssE,datay4ssEstdDev./2,'.black');
    
    gcoefsSine1 = sinFit(datay1ssE,datax1ssE); 
    gcoefsSine2 = sinFit(datay2ssE,datax2ssE); 
    gcoefsSine3 = sinFit(datay3ssE,datax3ssE); 
    gcoefsSine4 = sinFit(datay4ssE,datax4ssE); 
    
    
    fgsineDamp = @(p,x)(p(1).*exp(-p(2).*x).*sin(p(3).*x+p(4))+p(5));
    
    [gcoefsSineDamp1, gcoefsSineDampError1] = sinExpDampFit(datay1ssE,datax1ssE); 
    [gcoefsSineDamp2, gcoefsSineDampError2] = sinExpDampFit(datay2ssE,datax2ssE); 
    [gcoefsSineDamp3, gcoefsSineDampError3] = sinExpDampFit(datay3ssE,datax3ssE); 
    [gcoefsSineDamp4, gcoefsSineDampError4] = sinExpDampFit(datay4ssE,datax4ssE); 
    
    %frequencies:
    freq1 = gcoefsSineDamp1(3)/(10^(-3))/(2*pi);
    freq1min = gcoefsSineDampError1(3,1)/(10^(-3))/(2*pi);
    freq1max = gcoefsSineDampError1(3,2)/(10^(-3))/(2*pi);
    freq2 = gcoefsSineDamp2(3)/(10^(-3))/(2*pi);
    freq2min = gcoefsSineDampError2(3,1)/(10^(-3))/(2*pi);
    freq2max = gcoefsSineDampError2(3,2)/(10^(-3))/(2*pi);
    freq3 = gcoefsSineDamp3(3)/(10^(-3))/(2*pi);
    freq3min = gcoefsSineDampError3(3,1)/(10^(-3))/(2*pi);
    freq3max = gcoefsSineDampError3(3,2)/(10^(-3))/(2*pi);
    freq4 = gcoefsSineDamp4(3)/(10^(-3))/(2*pi);
    freq4min = gcoefsSineDampError4(3,1)/(10^(-3))/(2*pi);
    freq4max = gcoefsSineDampError4(3,2)/(10^(-3))/(2*pi);
    
    plot(fgsineDamp(gcoefsSineDamp1,1:100));
    plot(fgsineDamp(gcoefsSineDamp2,1:100),'r');
    plot(fgsineDamp(gcoefsSineDamp3,1:100),'g');
    plot(fgsineDamp(gcoefsSineDamp4,1:100),'black');
    
    plot(fgsine(gcoefsSine1,1:100));
    plot(fgsine(gcoefsSine2,1:100),'r');
    plot(fgsine(gcoefsSine3,1:100),'g');
    plot(fgsine(gcoefsSine4,1:100),'black');
    plot(datasx1Sorted,datasy1Sorted,'.');
    hold on;
    %plot(fgsine(gcoefsSine1,sort(datasx1)));
    plot(varData(71:140,1),sigmaR(71:140),'.r');
    plot(varData(141:210,1),sigmaR(141:210),'.g');
    plot(varData(211:end,1),sigmaR(211:end),'.black');
    
     plot(varData(set1,1),sigmaX(set1),'.r'); hold on; plot(varData(set1,1),sigmaY(set1),'.');
end

%------------ PCA Code ------------%

if(0)
        imagesToReshape = [];
        imagesToReshape = imageArrayC(:,:,1:35); %single images
        imagesToReshape = imagesAvgSet1;
        %imagesToReshape = imagesAvgSet2;
        %imagesToReshape = imagesAvgSet3;
        %imagesToReshape = imagesAvgSet4;
        
        numberOfImages = length(imagesToReshape(1,1,:));
        vectorLength = length(imagesToReshape(:,1,1))*length(imagesToReshape(1,:,1));
        
        %reshape images into vectors:
        vectorImages = []; vectorSum = zeros(vectorLength,1);
        for i=1:length(imagesToReshape(1,1,:))
            vectorImages(:,i) = reshape(imagesToReshape(:,:,i),[],1);
            vectorSum = vectorSum + vectorImages(:,i);
        end
        
        %Mean Image:
        meanImageVector = (1/numberOfImages).*vectorSum;
        undividedImage = reshape(vectorSum,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        meanImage = reshape(meanImageVector,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        
        normImages = [];
        for i=1:length(imagesToReshape(1,1,:))
            normImages(:,i) = vectorImages(:,i) - meanImageVector;
        end
        
        %N x p Bmatrix (N = number of images, p = pixels per image)
        Bmatrix = normImages;
        
        %[coeff, score, latent] = pca(ingredients)
        %each column of score corresponds to one principal component
        %latent stores the variances of the N principal components
        [Y01,P01,E01,tsquared,percentV] = pca(Bmatrix);
              
        
        %Reconstruct images using desired principle components
        eigenVectorBmatrix = []; imagesFromEVector = [];
        for i=1:numberOfImages
            eigenVectorBmatrix = P01(:,i);
            imagesFromEVector(:,:,i) = reshape(eigenVectorBmatrix,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        end
        
        
        %set of eigenmodes (principal components) large figure:
        figure(1);
        %set(gca,'XTickLabelMode', 'manual','XTickLabel', []);
        %set(gca,'YTickLabelMode', 'manual','YTickLabel', []);
        for i=1:numberOfImages
            subplot(floor(numberOfImages/5),5,i);
            imagesc(imagesFromEVector(:,:,i));
            set(gca,'XTickLabelMode', 'manual','XTickLabel', []);
            set(gca,'YTickLabelMode', 'manual','YTickLabel', []);         
        end
        
        %Fit to find frequency:
        %Mode 1:
        [gcoefsPCAmode1,gcoefsPCAerror1] = sinExpDampFit(Y01(:,1),datax1ssE);
        freqPCAmode1 = gcoefsPCAmode1(3)/(10^(-3))/(2*pi);
        freqPCAmode1ErrorMin = gcoefsPCAerror1(3,1)/(10^(-3))/(2*pi);
        freqPCAmode1ErrorMax = gcoefsPCAerror1(3,2)/(10^(-3))/(2*pi);
        
        figure(7);
        subplot(2,2,1);
        imagesc(imagesFromEVector(:,:,1));
        subplot(2,2,2);
        plot(fgsineDamp(gcoefsPCAmode1,1:100)); hold on;
        plot(datax1ssE,Y01(:,1),'.'); 
        text(40,0.32,['\omega = 2\pi \times ' num2str(freqPCAmode1,4) ' (' num2str(freqPCAmode1ErrorMin,4) ',' num2str(freqPCAmode1ErrorMax,4) ')']); 
        hold off;
        
        %Mode 2:
        [gcoefsPCAmode2,gcoefsPCAerror2] = sinExpDampFit(Y01(:,2),datax1ssE);
        freqPCAmode2 = gcoefsPCAmode2(3)/(10^(-3))/(2*pi);
        freqPCAmode2ErrorMin = gcoefsPCAerror2(3,1)/(10^(-3))/(2*pi);
        freqPCAmode2ErrorMax = gcoefsPCAerror2(3,2)/(10^(-3))/(2*pi);
 
        %figure(3);
        subplot(2,2,3);
        imagesc(imagesFromEVector(:,:,2));
        subplot(2,2,4);        
        plot(fgsineDamp(gcoefsPCAmode2,1:100)); hold on;
        plot(datax1ssE,Y01(:,2),'.');
        text(40,0.32,['\omega = 2\pi \times ' num2str(freqPCAmode2,4) ' (' num2str(freqPCAmode2ErrorMin,4) ',' num2str(freqPCAmode2ErrorMax,4) ')']); 
        hold off;
        
        %Mode 3:
        [gcoefsPCAmode3,gcoefsPCAerror3] = sinExpDampFit(Y01(:,3),datax1ssE);
        freqPCAmode3 = gcoefsPCAmode3(3)/(10^(-3))/(2*pi);
        freqPCAmode3ErrorMin = gcoefsPCAerror3(3,1)/(10^(-3))/(2*pi);
        freqPCAmode3ErrorMax = gcoefsPCAerror3(3,2)/(10^(-3))/(2*pi);

        subplot(3,2,5);
        imagesc(imagesFromEVector(:,:,3));
        subplot(3,2,6);
        plot(fgsineDamp(gcoefsPCAmode3,1:100)); hold on;
        plot(datax1ssE,Y01(:,3),'.');
        text(40,0.32,['\omega = 2\pi \times ' num2str(freqPCAmode3,4) ' (' num2str(freqPCAmode3ErrorMin,4) ',' num2str(freqPCAmode3ErrorMax,4) ')']);
        hold off;
        
end


%----- Values ------%

if(0)

    atomNumsRun = [45000, 20000, 15000, 10000];
    %Radial Fitting Results%
    breathingModeFreqsGaussFit = [freq1, freq2, freq3, freq4];
    breathingModeFreqsGaussFitMin = [freq1min, freq2min, freq3min, freq4min];
    breathingModeFreqsGaussFitMax = [freq1max, freq2max, freq3max, freq4max];
    %PCA Results% 
    dipoleModeFreqs = [20.7, 21.83, 20.42, 21.18];
    dipoleModeFreqsMin = [19.98, 20.93, 19.22, 19.66];
    dipoleModeFreqsMax = [21.4, 22.73, 21.61, 22.7];
    breathingModeFreqs = [38.23, 38.78, 39.15, 42.81];
    breathingModeFreqsMin = [37.99, 38.37, 38.56, 40.98];
    breathingModeFreqsMax = [38.48, 39.18, 39.75, 44.65];
    
    figure(10); 
    errorbar(atomNumsRun,dipoleModeFreqs,dipoleModeFreqs-dipoleModeFreqsMin,dipoleModeFreqsMax-dipoleModeFreqs,'.');  
    
    figure(11);
    errorbar(atomNumsRun,breathingModeFreqsGaussFit,breathingModeFreqsGaussFit-breathingModeFreqsGaussFitMin,breathingModeFreqsGaussFitMax-breathingModeFreqsGaussFit,'.r');
    hold on;
    errorbar(atomNumsRun,breathingModeFreqs,breathingModeFreqs-breathingModeFreqsMin,breathingModeFreqsMax-breathingModeFreqs,'.');
    hold off;

end


if(0)
    %plot strings:
    plot(pixelNumbers,widthsR,'.');
    plot(motFets,widthsR,'.');
    plot(varData(:,1),sigmaR,'.');
    plot(pixelNumbers,widthsX,'.'); hold on; plot(pixelNumbers,widthsY,'.r');
    plot(nROIlogfile,widthXlogfile,'.'); hold on; plot(nROIlogfile,widthYlogfile,'.r');
    plot(nROIlogfile,widthXlogfile,'.r'); hold on; plot(pixelNumbers,widthsX,'.');
    plot(nROIlogfile,widthXlogfile,'.r'); hold on; plot(pixelCounts,sigmaX,'.');
    plot(pixelCounts,'.'); hold on; plot(nROIlogfile','.r');
end