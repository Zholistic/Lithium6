directory = 'C:\Data\150903_Top_Side_Compare_Magnification\150903_magnification_compare_top\';
directory = 'C:\Data\150904_Top_Side_Compare_Magnification_600usTOF\150904_magnification_compare_top\';
directory = 'C:\Data\150904_Top_Side_Compare_Magnification_600usTOF\150904_5ForCalib\';
directory = 'C:\Data\150906_mag_side\';
date = '150903';
date = '150904';
date = '150906';
camera = 'sidecam';
varstring = 'imagenumber';
magfield = '832p2G';
bins = 36;
%varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %2.84 um topcam, topcam magnification = 4.58
pixelLength = 3.75e-6 / 1.9; %1.9 magnification (old 1.4 mag) ??
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
Isat = 10^6;
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
    imageArray(:,:,i) = PullFTS(fileLocList{i},raw,Isat);
end

%Crop images:
%Crop images:
imageArrayC = []; imageArrayTC = [];
ROIx = 200:980;
ROIy = 500:675;
%The cross is inside the region specified above. 
%Can improve the error by tightening the Cross ROI up
%for that particular dataset. 
CrossROIx = 271:498; 
CrossROIy = 73:85;
TightROIy = 565:580;
TightROIx = 538:619;

imageArrayC = imageArray(ROIy,ROIx,:);
imageArrayTC = imageArray(TightROIy,TightROIx,:);


%Display every X image:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayC(:,:,i));        
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayTC(:,:,i));        
    end
end
end

%%%%%Fit Gaussians:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = []; sigmaX = [];
sigmaY = [];
for i=1:length(imageArrayC(1,1,:))
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    
    %Shift wings to zero:
    shiftFactor = (gcoefsXi(4,i)+gcoefsYi(4,i))/2;
    imageArrayC(:,:,i) = imageArrayC(:,:,i) - shiftFactor;
    imageArrayTC(:,:,i) = imageArrayTC(:,:,i) - shiftFactor;
    
    %Refit:
    gcoefsX(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    gcoefsY(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))
    
    centers(:,i) = [ceil(gcoefsX(2,i)), ceil(gcoefsY(2,i))]; %center = [x y]
    sigmaX(:,i) = gcoefsX(3,i);
    sigmaY(:,i) = gcoefsY(3,i);
end


%%%%%Second moment:
xvector = [];
xvector = 1:length(imageArrayTC(1,:,1));
COMx = [];
for i=1:length(imageArrayTC(1,1,:))
        %Sum of Intensity*pixel location / sum of intensity
        COMx(i) = sum(mean(imageArrayTC(:,:,i),1).*xvector) / sum(mean(imageArrayTC(:,:,i),1));
end

%Second moment x direction
SMomX = []; sigmaXsm = [];
for i=1:length(imageArrayTC(1,1,:))
    SMomX(i) = sum(mean(imageArrayTC(:,:,i),1).*(xvector - COMx(i)).^2) / sum(mean(imageArrayTC(:,:,i),1));
end

sigmaXsm = sqrt(SMomX);

yvector = [];
yvector = 1:length(imageArrayTC(:,1,1));
COMy = [];
for i=1:length(imageArrayTC(1,1,:))
        %Sum of Intensity*pixel location / sum of intensity
        COMy(i) = sum(mean(imageArrayTC(:,:,i),2)'.*yvector) / sum(mean(imageArrayTC(:,:,i),2));
end

%Second moment y direction
SMomY = []; sigmaYsm = [];
for i=1:length(imageArrayTC(1,1,:))
    SMomY(i) = sum(mean(imageArrayTC(:,:,i),2)'.*(yvector - COMy(i)).^2) / sum(mean(imageArrayTC(:,:,i),2));
end

sigmaYsm = sqrt(SMomY);


%%%%%Temperatures:
if(0)
TonTFs = [];
for i=1:length(imageArrayC(1,1,:))
    TonTFs(i) = 1/(log(1+exp(gcoefsPolyLog1(2,i)/gcoefsPolyLog1(3,i)^2)));
end
end


%Display every X fit:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        plot(fglz(gcoefsX(:,i),1:400));hold on; plot(mean(imageArrayC(CrossROIy,:,i),1),'r'); hold off;      
    end
end
%These ARE the plots you're looking for:
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        plot(fglz(gcoefsY(:,i),1:180));hold on; plot(mean(imageArrayC(:,CrossROIx,i),2),'r'); hold off;      
    end
end
end

%%%%%Atom numbers:
%Tight ROI array:
pixelCounts = [];
for i=1:length(imageArrayTC(1,1,:))
    pixelCounts(i) = sum(sum(imageArrayTC(:,:,i)));
end

%Pixel to Atom Correction:
pixelCounts(:) = varData(:,5)*0.42;


%Sort varData:
indexs = [];
[sortedVarData,indexs] = sort(varData);
%indexs(:,1) is a vector of the sort.

pixelCountsSort = []; sigmaXSort = []; sigmaYSort = [];
TonTFsSort = []; imageArrayCSort = [];
sigmaRsmSort = []; TonTFsSort = [];

for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
    sigmaXsMomSort(i) = sigmaXsm(indexs(i));
    sigmaYsMomSort(i) = sigmaYsm(indexs(i));
    
    imageArrayCSort(:,:,i) = imageArrayC(:,:,indexs(i));
end

%Average over same BField data points:
j=1; runTotal = 0; magFields = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
imageArrayAvgs = [];
sMomentXStdDev = []; sMomentYStdDev = []; prev = sortedVarData(1,1);
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
        magFields(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        sMomentX(j) = mean(sigmaXsMomSort(i-runTotal:i));
        sMomentY(j) = mean(sigmaYsMomSort(i-runTotal:i));
        sMomentXStdDev(j) = std(sigmaXsMomSort(i-runTotal:i));
        sMomentYStdDev(j) = std(sigmaYsMomSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        runTotal = 0;
        j = j+1;
    end
    if( i == length(sortedVarData))
        %last run
        %disp('hit last run')
        if(i == runTotal)
            runTotal = runTotal-1;
        end
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        magFields(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        sMomentX(j) = mean(sigmaXsMomSort(i-runTotal:i));
        sMomentY(j) = mean(sigmaYsMomSort(i-runTotal:i));
        sMomentXStdDev(j) = std(sigmaXsMomSort(i-runTotal:i));
        sMomentYStdDev(j) = std(sigmaYsMomSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        runTotal = 0;
        j = j+1;
    end
    
    prev = curr;
end
%Fit functions to the averaged images:
gcoefsXa = []; gcoefsYa = []; gcoefsYaError = []; gcoefsXaError = [];
for i=1:length(imageArrayAvgs(1,1,:))
    [gcoefsXa(:,i),gcoefsXaError(:,:,i)] = gausFit1DLockZero(mean(imageArrayAvgs(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    [gcoefsYa(:,i),gcoefsYaError(:,:,i)] = gausFit1DLockZero(mean(imageArrayAvgs(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))       
end

%function:
fg1d = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2))); 

widthsYavg = gcoefsYa(3,:);
widthsXavg = gcoefsXa(3,:);

widthsYavgError = []; widthsXavgError = [];
for i=1:length(widthsYavg)
    widthsYavgError(i) = widthsYavg(i) - gcoefsYaError(3,1,i);
    widthsXavgError(i) = widthsXavg(i) - gcoefsXaError(3,1,i);
end

%Images:
if(0)
for i=1:length(imageArrayAvgs(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        plot(fg1d(gcoefsXa(:,i),1:800));hold on; plot(mean(imageArrayAvgs(CrossROIy,:,i),1),'r'); hold off;      
    end
end
for i=1:length(imageArrayAvgs(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayAvgs(:,:,i));        
    end
end
%Raw Images
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayC(:,:,i));        
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayC(CrossROIy,CrossROIx,i));        
    end
end
%X Profiles:
for i=1:10
    figure(i);
    plot(mean(imageArrayC(CrossROIy,:,i),1))
end
end





