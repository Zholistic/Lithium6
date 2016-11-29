directory = 'C:\Data\EOS_Data\140908_2d_eos_880G_10k_atoms_750ms_laser_ramp_0p25_10isat\';
%directory = 'C:\Data\EOS_Data\140910_2d_eos_972G_10us_0p25_isat_1us_10isat_750ms_laser_ramp\';
date = '140908';
camera = 'top';
varstring = 'isat';
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

imageArray = []; imageState1Array = []; beamImageArray = [];
%Pull images:
for i=1:length(fileLocList)
    Isat = varData(i,1);
    Isat = 134;
    if(Isat < 1) 
        Isat = 1340;
    end
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArray(:,:,i) = atom2Image(:,:,1);
    imageState1Array(:,:,i) = atom1Image(:,:,1);
    beamImageArray(:,:,i) = beamImage(:,:,1);    
end

%Crop images:
imageArrayC = []; imageArrayTC = [];
imageArrayHighInt = []; imageArrayLowInt = [];
ROIx = 1:length(imageArray(1,:,1));
ROIy = 1:length(imageArray(:,1,1));
CrossROIy = 30:160; 
CrossROIx = 40:180; %The cross is inside the region specified above.
CrossROIy = 1:100;
CrossROIx = 1:110;

TightROIx = 50:170;
TightROIy = 30:140;
%Split into high and low intensity arrays
imageArrayC = imageArray(TightROIy,TightROIx,:);
imageArrayTC = imageArray(TightROIy,TightROIx,:);

imageArrayHighInt = imageArrayC(:,:,1:200);
imageArrayLowInt = imageArrayC(:,:,201:end);

%Display every X image:
if(0)
for i=1:length(imageArray(1,1,:))
    if(mod(i,10) == 0)       
        figure(i);
        imagesc(imageArray(:,:,i));        
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,10) == 0)       
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
    radProfiles(:,:,i) = radProfilesT(:,1:end-5,i);
end

%Shift wings of radial profiles to zero:
shiftBy = [];
for i=1:length(radProfiles(1,1,:))
    shiftBy(i) = mean(radProfiles(1,72:84,i));
    radProfiles(1,:,i) = radProfiles(1,:,i) - shiftBy(i);
end

if(0)
    for i=1:length(imageArrayC(1,1,:))
        if(mod(i,2) == 0)
            figure(i);
            plot(fgr(gcoefsR(:,i),1:180),'g'); hold on;  plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); hold off;
        end
    end
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
    %imageArrayC(:,:,i) = imageArrayC(:,:,i) - shiftFactor(i);
        
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
    plot(varData(1:30,1),sigmaX(1:30),'.'); hold on; plot(varData(31:end,1),sigmaX(31:end),'.r');
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


%Re-order each array:
sortedVarData = []; indexs = []; sigmaRSort = [];
sigmaXSort = []; sigmaYSort = []; sigmaR2p5Sort = [];
sigmaRsmSort = []; TonTFsSort = []; imageArrayCSort = [];

reorder = 0;
if(reorder)
    [sortedVarData,indexs] = sort(varData);
%indexs(:,1) is a vector of the sort (varstring).
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
else
    sortedVarData = varData;
    for i=1:length(sigmaX)
        sigmaXSort(i) = sigmaX(i);
        sigmaYSort(i) = sigmaY(i);
        pixelCountsSort(i) = pixelCounts(i);
        sigmaRSort(i) = sigmaR(i);
        %sigmaR2p5Sort(i) = sigmaR2p5(indexs(i));
        %sigmaRsmSort(i) = sigmaRsm(indexs(i));
        TonTFsSort(i) = TonTFs(i);
        imageArrayCSort(:,:,i) = imageArrayC(:,:,i);
    end
end


%Average over same motfet data points:
%Bug fix history: 150109 fixed endNum
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
widthsR = []; stdDevWidthsR = []; widthsR2p5 = [];
widthsPsm = []; stdDevWidthsPsm = []; TonTFsm = []; stdDevTonTFsm = [];
sMomentXStdDev = []; sMomentYStdDev = []; stdDevWidthsR2p5 = [];
imageArrayAvgsNoCenter = [];
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
        imageArrayAvgsNoCenter(:,:,j) = mean(imageArrayCSort(:,:,i-runTotal:endNum),3);
        
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
        imageArrayAvgsNoCenter(:,:,j) = mean(imageArrayCSort(:,:,i-runTotal:i),3);
        
        motFets(j) = sortedVarData(i);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        runTotal = 1;
        j = j+1;
    end
    
    prev = curr;
end

%Radially averaged average profiles:
radProfilesAvgReorder = []; radProfilesTAvg = []; centerAvgs = [];
disp('Radially averaging...');
for i=1:length(imageArrayAvgs(1,1,:))
    [radProfilesTAvg(:,:,i),centerAvgs(:,i)] = radAverageReorder(imageArrayAvgs(:,:,i));
    radProfilesAvgReorder(:,:,i) = radProfilesTAvg(:,1:end-5,i);
end

radProfilesAvg = [];
n=170;
radProfilesAvg(2,:,1) = meanNelements(radProfilesAvgReorder(2,:,1),n);
radProfilesAvg(1,:,1) = meanNelements(radProfilesAvgReorder(1,:,1),n);
radProfilesAvg(2,:,2) = meanNelements(radProfilesAvgReorder(2,:,2),n);
radProfilesAvg(1,:,2) = meanNelements(radProfilesAvgReorder(1,:,2),n);

%EOSGenerate(radProfilesAvg(1,:,1)./(kpixelLength^2),radProfilesAvg(2,:,1),2:65,omegaRVector(2))

atomNumBC = [];
atomNumBC(1) = sum(sum(imageArrayAvgs(:,:,1)));
atomNumBC(2) = sum(sum(imageArrayAvgs(:,:,2)));

pixel2NROI = 1/(atomNumBC(1)/mean(varData(1:200,5)));
NROI2HighInt = atomNumBC(2)/(atomNumBC(1)*pixel2NROI);
HighInt2TwoSpinState = 2;

atomCorr = pixel2NROI*NROI2HighInt*HighInt2TwoSpinState;
correctedNum = atomNumBC(1)*pixel2NROI*NROI2HighInt*HighInt2TwoSpinState;

radProfilesAvg(1,:,1) = radProfilesAvg(1,:,1).*atomCorr;

%Bulk:

%Radially averaged profiles:
radProfilesToReorder = []; radProfilesTT = []; center1 = []; radProfilesIndvToEOS = [];
disp('Radially averaging...');
for i=1:length(imageArrayC(1,1,:))
    [radProfilesTT(:,:,i),center1(:,i)] = radAverageReorder(imageArrayC(:,:,i));
    radProfilesToReorder(:,:,i) = radProfilesTT(:,1:end-5,i);
    radProfilesIndvToEOS(1,:,i) = meanNelements(radProfilesToReorder(1,:,i),300);
    radProfilesIndvToEOS(2,:,i) = meanNelements(radProfilesToReorder(2,:,i),300);
end

radProfilesIndvToEOS(1,:,1:200) = radProfilesIndvToEOS(1,:,1:200).*atomCorr;

%Shift wings of radial profiles to zero:
%shiftBy = [];
%for i=1:length(radProfilesIndvToEOS(1,1,:))
    %shiftBy(i) = mean(radProfilesIndvToEOS(1,39:66,i));
    %radProfilesIndvToEOS(1,:,i) = radProfilesIndvToEOS(1,:,i) - shiftBy(i);
    %radProfilesIndvToEOS(1,:,i) = radProfilesIndvToEOS(1,:,i) - min(radProfilesIndvToEOS(1,:,i));
%end

if(0)
    for i=1:length(radProfilesIndvToEOS(1,1,:))
        if(mod(i,10) == 0)
            figure(i);
            plot(radProfilesIndvToEOS(2,:,i),radProfilesIndvToEOS(1,:,i),'r'); hold off;
        end
    end
end

calcRegion = 2:45;
a0 = 5.29e-11; %o.0
a2dVector = a0.*[47194.1 66130.5 141838 302318];
a2dVectorFiveOneFive = a0.*[52679 75445.1 169680 377935];
a2d = a2dVectorFiveOneFive(2);
omegaRVector = 2.*pi.*[24.88224858 25.09885405 25.66767827 26.3891196];
omegaR = omegaRVector(2);
saveOn = 0;
smoothOn = 1;
zeroOn = 0;
savestring = 'indv clouds smoothed';
savestring = 'test_disregard';
%EOSGenerateBulk(radProfilesIndvToEOS(1,:,1:200)./(kpixelLength^2),radProfilesIndvToEOS(2,:,1:200),calcRegion,omegaR,880,a2d)
%EOSGenerateBulk(radProfilesIndvToEOS(1,:,1:200)./(kpixelLength^2),radProfilesIndvToEOS(2,:,1:200),calcRegion,omegaR,880,a2d,smoothOn,zeroOn,saveOn,savestring)
%single in bulk function:
%EOSGenerateBulk(radProfilesAvg(1,:,1)./(kpixelLength^2),radProfilesAvg(2,:,1),2:65,omegaR,880,a2d,smoothOn,zeroOn,saveOn, savestring);

calcRegion = 2:65;
omegaZ = 5150*2*pi;
cutoffLeft = 23;
cutoffRight = 0;
%fitVirial(radProfilesAvg(1,:,1)./(kpixelLength^2), radProfilesAvg(2,:,1),calcRegion, omegaZ, omegaR, 880, a2d, smoothOn, zeroOn, cutoffLeft, cutoffRight)

%Shift wings of radial profiles to zero:
%radProfilesAvg(1,:,1) = radProfilesAvg(1,:,1) - min(radProfilesAvg(1,:,1));
shiftByAvgs = [];
for i=1:length(radProfilesAvg(1,1,:))
    shiftByAvgs(i) = mean(radProfilesAvg(1,51:61,i));
    radProfilesAvg(1,:,i) = radProfilesAvg(1,:,i) - shiftByAvgs(i);
end
%radProfilesAvg(1,:,i) = radProfilesAvg(1,:,i) - min(radProfilesAvg(1,:,1)); !!!!


%Shift images by found factor:
imageArrayAvgsZerod = [];
for i=1:length(imageArrayAvgs(1,1,:))
    imageArrayAvgsZerod(:,:,i) = imageArrayAvgs(:,:,i) - shiftByAvgs(i);
end


if(0)
    for i=1:length(imageArrayAvgs(1,1,:))
        if(mod(i,1) == 0)
            figure(i);
            hold on; plot(radProfilesAvg(2,:,i),radProfilesAvg(1,:,i),'r'); line([0 100],[0 0]); hold off;
        end
    end
end

%%%%%%%%%%%%%%%%%----- Radially averaged density vs potential (individual)
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
radProfilesToParse = [];
omegaRVector = 2.*pi.*[24.88224858 25.09885405 25.66767827 26.3891196];
magVector = [865 880 920 972];
kpixelLength = (13e-6*(83/400)); %kristian pixel length


radProfilesToParse = radProfiles(:,:,1:200);

%Zeroing
shiftByAvgs = [];
for i=1:length(radProfilesToParse(1,1,:))
    shiftByAvgs(i) = mean(radProfilesToParse(1,59:70,i));
    radProfilesToParse(1,:,i) = radProfilesToParse(1,:,i) - shiftByAvgs(i);
end

if(0)
    for i=1:length(radProfilesLowInt(1,1,:))
        if(mod(i,10) == 0)
            figure(i);
            hold on; line([0 100],[0 0]); plot(radProfilesToParse(2,:,i),radProfilesToParse(1,:,i),'r'); hold off;
        end
    end
end

%prepare the x-axis (potential in nk)
radiusVectors = []; radiusVectorsMeters = []; potentials = []; potentials_nk = [];
for i = 1:length(radProfilesToParse(2,1,:))
    radiusVectors(:,i) = 1:(length(radProfilesToParse(2,:,i)));
    radiusVectorsMeters(:,i) = radiusVectors(:,i).*kpixelLength;
    
    potentials(:,i) = 0.5 .* massL6 .* (omegaRVector(2)^2).* radiusVectorsMeters(:,i).^2; %Omega R to corresponding B Field
    potentials_nk(:,i) = (potentials(:,i)./kB).*(10^9);
end



%prepare the y-axis (number density)
radProfilesDensityAdjusted = []; radproftoFitReal = [];

for i=1:length(radProfilesToParse(1,1,:))
    radproftoFitReal(:,i) = radProfilesToParse(1,:,i)./(kpixelLength^2); %y values real units
    
    radProfilesDensityAdjusted(:,i) = radproftoFitReal(:,i).*2.*1.2; %2 for spin states and 1.2 for correction factor
end

%plot(potential_nk,radproftoFitReal);

if(0)
    for i=1:length(radProfilesToParse(1,1,:))
        if(mod(i,10) == 0)
            figure(i);
            subplot(1,2,1);
            imagesc(imageArrayC(:,:,i));
            subplot(1,2,2);
            hold on; line([0 600],[0 0]); plot(potentials_nk(:,i),radproftoFitReal(:,i),'.r'); hold off;
        end
    end
end

%kappa and p calculations:
calcRegion = 5:55;
kappaTilde = []; pTilde = []; p = [];
for i=1:length(radProfilesToParse(1,1,:))
    
    kappaTilde(:,i) = (-1)*(pi * hbar^2 / massL6) * gradient(radproftoFitReal(calcRegion,i),potentials(calcRegion,i));
           
    
    for j=1:length(calcRegion)
        p(j,i) = trapz(potentials(j:calcRegion(end),i),radproftoFitReal(j:calcRegion(end),i)); %integral from end to j point
        n = radproftoFitReal(j,i);     %n, density on j'th point (varies across the potential)
        pTilde(j,i) = ((2*massL6)/(n^2 * hbar^2 * pi)) * p(j,i);
    end

end


figure(1000);
hold on;
for i=1:length(pTilde(1,:))
    hLine = plot(pTilde(5:40,i),kappaTilde(5:40,i),'.');
    %scatter(pTilde(5:40,i),kappaTilde(5:40,i),'filled', ...
        %'MarkerFaceAlpha','MarkerFaceColor','blue')
    %plot(pTilde(5:50,i),'.');
    %hMarkers = hLine.MarkerHandle;
    %hMarkers.FaceColorData = uint8(255*[1;0;0;0.3]);  % Alpha=0.3 => 70% transparent red
end
axis([-4 2 -6 6])
%axis([0 50 -1.5 5])
hold off;

cloudNumber = 2;
figure(9); plot(pTilde(:,cloudNumber),kappaTilde(:,cloudNumber),'.');
title('kappa vs p');
figure(10); h1 = plot(pTilde(:,cloudNumber),'.');
title('p vs r (pixel)');
figure(11); plot(kappaTilde(:,cloudNumber),'.');
title('kappa vs r (pixel)');
figure(12); plot(potentials(calcRegion,cloudNumber),radproftoFitReal(calcRegion,cloudNumber),'.');
title('density vs potential');

%%%%%%%%%%%%%%%%%%%%----Fit Virial:
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
magfield = 865;
magfield = 972;
omegaz = 5789 * 2 * pi;
omegar = 25 * 2 * pi; %Double check this weak trapping

omegaRVector = 2.*pi.*[24.88224858 25.09885405 25.66767827 26.3891196];
a0 = 5.29e-11; %o.0
a2dVector = a0.*[47194.1 66130.5 141838 302318];
magVector = [865 880 920 972];
az = sqrt(hbar/(massL6*omegaz));

pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
kpixelLength = (13e-6*(83/400)); %kristian pixel length

radproftoFityReal = radProfilesAvg(1,:,1)./(kpixelLength^2);
%radproftoFityRealConv = radproftoFityReal.*10^(-11);
%radproftoFitxReal = radProfilesAvg(2,58:76,2).*pixelLength;
radiusVector = 1:(length(radProfilesAvg(2,:,1)));
radiusVectorMeters = radiusVector.*kpixelLength;
potential = 0.5 .* massL6 .* omegaRVector(2)^2 .* radiusVectorMeters.^2;
potential_nk = (potential./kB).*(10^9);
%plot(potential_nk,radProfilesAvg(1,:,i));

%convert radProfile pixels to potential:
radProfilesPotential = 0.5 .* massL6 .* omegaRVector(2)^2 .*(radProfilesAvg(2,:,1).*pixelLength).^2;
%plot(radProfilesPotential,radProfilesAvg(1,:,1)./(pixelLength^2));
%radProfilesPotentialConv = radProfilesPotential.*10^30;
radProfilesDensityAdjusted = radproftoFityReal.*2.*1.2; %2 for spin states and 1.2 for correction factor
%plot(potential_nk,radproftoFityReal);
%plot(potential,radproftoFityReal);

a0 = 5.29e-11; %o.0
a2d = 47186.7 * a0; %From mathematica spreadsheet TODO generate text file
a2d = 301946 * a0; %972G, 5.789*2pi hz
eb = (hbar^2)/(a2d^2 * massL6);


%Scale x&y to H.O. dimensionless 
%Density SI = m^-2
radProfilesDensityConv = radProfilesDensityAdjusted./(massL6*omegaz / hbar);
%Potential SI = kg * s^-2 * m^2 // Energy units = kg * m^2 * s^-2
radProfilesPotentialConv = radProfilesPotential./((hbar*omegaz)/2); 
% plot(radProfilesPotentialConv,radProfilesDensityConv);

%kappa and p calculations:
calcRegion = 2:60;
kappaTildeAvg = []; pTildeAvg = []; pAvg = []; dydx = [];

%kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) * gradient(radProfilesDensityAdjusted(calcRegion),potential(calcRegion));
dydx = diff([eps radProfilesDensityAdjusted(calcRegion)])./diff([eps potential(calcRegion)]);
%kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) * gradient(radProfilesDensityAdjusted(calcRegion),potential(calcRegion));
kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) .* dydx;


for j=1:length(calcRegion)
    pAvg(j) = trapz(potential(j:calcRegion(end)),radProfilesDensityAdjusted(j:calcRegion(end))); %integral from end to j point
    n = radProfilesDensityAdjusted(j);     %n, density on j'th point (varies across the potential)
    pTildeAvg(j) = ((2*massL6)/(n^2 * hbar^2 * pi)) * pAvg(j);
end

close all;
%averaged cloud kappa and p's
figure(8); plot(radiusVectorMeters,radProfilesDensityAdjusted,'.');
title('Density Adjusted vs radius (meters)');
figure(9); plot(pTildeAvg,kappaTildeAvg,'.');
title('kappaTilde vs pTilde'); axis([0 14 0 2]);
figure(10); h1 = plot(potential(calcRegion),pTildeAvg,'.');
title('pTilde vs potential (J)');
figure(11); plot(potential(calcRegion),kappaTildeAvg,'.');
title('kappaTilde vs potential (J)');
figure(12); plot(potential_nk(calcRegion),radProfilesDensityAdjusted(calcRegion),'.');
title('density vs potential');

fitrange = 54:64;

virialRsquares = []; tempGuess = []; virialTreturned = []; virialMu0returend = [];
virialFitFuncs = [];
for i=1:60;
tempGuess(i) = 10e-9 + i*(0.5e-9); %Guess 20nK
betaeb = eb/(kB*tempGuess(i));

b2 = huicoeffs(betaeb,1,0);
b3 = huicoeffs(betaeb,2,0);

[f1, gof, foutput] = virial2Fit(radProfilesDensityConv(fitrange),radProfilesPotentialConv(fitrange),[b2 b3]);
%virial2Fit(radProfilesDensityAdjusted(fitrange),radProfilesPotential(fitrange),[b2 b3]);

%figure(i);
%plot(f1); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.'); hold off;

virialRsquares(i) = getfield(gof,'rsquare');
virialTreturned(i) = getfield(f1,'a');
virialMu0returend(i) = getfield(f1,'b');
virialFitFuncs{i} = f1;
end

if(0)
    close all;
    plot(tempGuess,'.'); hold on; plot(virialTreturned,'.r');
    figure(2); plot(virialRsquares,'.');
    
    figure(3); plot(virialFitFuncs{19}); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.');
end

%T / TF calc:

n = 0.25*(massL6*omegaz / hbar)/2;
EFermi = (hbar^2 * 4 * pi * n) / (2*massL6);
TFermi = EFermi/kB;

TonTF = 20.24e-9 / TFermi;

atomNum = [];
for i=1:length(imageArrayAvgs(1,1,:))
    atomNum(i) = sum(sum(imageArrayAvgs(:,:,i)));
end