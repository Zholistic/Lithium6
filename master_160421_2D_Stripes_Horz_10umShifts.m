directory = 'C:\Data\160421_stripes2D_shift_horizontal\';
%directory = 'C:\Data\151003_New_breath_test\';
%directory = 'C:\Data\150904_magnification_check3_top\';
%NOT ACTUALLY 724G. like 780 or 750, check lab book...
date = '160421';
camera = 'topcam';
varstring = 'BraggFrequency';
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
    OD = 1;
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

minusTwenty = 1:10;
minusTen = 11:20;
zero = 21:30;
plusTen = 31:40;
plusTwenty = 41:50;

imageArray1 = imageArrayC(:,:,minusTwenty);
imageArray2 = imageArrayC(:,:,minusTen);
imageArray3 = imageArrayC(:,:,zero);
imageArray4 = imageArrayC(:,:,plusTen);
imageArray5 = imageArrayC(:,:,plusTwenty);

imageArray1Avg = mean(imageArray1(:,:,:),3);
imageArray2Avg = mean(imageArray2(:,:,:),3);
imageArray3Avg = mean(imageArray3(:,:,:),3);
imageArray4Avg = mean(imageArray4(:,:,:),3);
imageArray5Avg = mean(imageArray5(:,:,:),3);

fig80 = figure(80); imagesc(imageArray1Avg);
fig90 = figure(90); imagesc(imageArray2Avg);
fig100 = figure(100); imagesc(imageArray3Avg);
fig110 = figure(110); imagesc(imageArray4Avg);
fig120 = figure(120); imagesc(imageArray5Avg);

F(1) = getframe(fig80);
F(2) = getframe(fig90);
F(3) = getframe(fig100);
F(4) = getframe(fig110);
F(5) = getframe(fig120);

fig1 = figure(1); movie(fig1,F,20,5);

imageArray1AvgRot = imrotate(imageArray1Avg,-20);
imageArray2AvgRot = imrotate(imageArray2Avg,-20);
imageArray3AvgRot = imrotate(imageArray3Avg,-20);
imageArray4AvgRot = imrotate(imageArray4Avg,-20);
imageArray5AvgRot = imrotate(imageArray5Avg,-20);

fig80 = figure(80); imagesc(imageArray1AvgRot);
fig90 = figure(90); imagesc(imageArray2AvgRot);
fig100 = figure(100); imagesc(imageArray3AvgRot);
fig110 = figure(110); imagesc(imageArray4AvgRot);
fig120 = figure(120); imagesc(imageArray5AvgRot);

xregionToAvg = 125:151;
yregionToAvg = 38:230;

image1Slice = mean(imageArray1AvgRot(yregionToAvg,xregionToAvg),2);
image2Slice = mean(imageArray2AvgRot(yregionToAvg,xregionToAvg),2);
image3Slice = mean(imageArray3AvgRot(yregionToAvg,xregionToAvg),2);
image4Slice = mean(imageArray4AvgRot(yregionToAvg,xregionToAvg),2);
image5Slice = mean(imageArray5AvgRot(yregionToAvg,xregionToAvg),2);


figure(1); subplot(5,1,1); plot(image1Slice);
subplot(5,1,2); plot(image2Slice);
subplot(5,1,3); plot(image3Slice);
subplot(5,1,4); plot(image4Slice);
subplot(5,1,5); plot(image5Slice);

gcoefsI = [];
gcoefsI(:,1) = gausFit1D(image1Slice);
gcoefsI(:,2) = gausFit1D(image2Slice);
gcoefsI(:,3) = gausFit1D(image3Slice);
gcoefsI(:,4) = gausFit1D(image4Slice);
gcoefsI(:,5) = gausFit1D(image5Slice);

%zeroing
image1Slice = image1Slice - mean(image1Slice(1:18));
image2Slice = image2Slice - mean(image2Slice(1:18));
image3Slice = image3Slice - mean(image3Slice(1:18));
image4Slice = image4Slice - mean(image4Slice(1:18));
image5Slice = image5Slice - mean(image5Slice(1:18));

gcoefs = [];
gcoefs(:,1) = gausFit1D(image1Slice);
gcoefs(:,2) = gausFit1D(image2Slice);
gcoefs(:,3) = gausFit1D(image3Slice);
gcoefs(:,4) = gausFit1D(image4Slice);
gcoefs(:,5) = gausFit1D(image5Slice);


fg1d = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));

xs = 1:length(image1Slice);

figure(1); plot(image1Slice); hold on; plot(fg1d(gcoefs(:,1),xs),'r');
figure(2); plot(image2Slice); hold on; plot(fg1d(gcoefs(:,2),xs),'r');
figure(3); plot(image3Slice); hold on; plot(fg1d(gcoefs(:,3),xs),'r');
figure(4); plot(image4Slice); hold on; plot(fg1d(gcoefs(:,4),xs),'r');
figure(5); plot(image5Slice); hold on; plot(fg1d(gcoefs(:,5),xs),'r');

figure(10); plot(image1Slice' - fg1d(gcoefs(:,1),xs));
figure(11); plot(image2Slice' - fg1d(gcoefs(:,2),xs));
figure(12); plot(image3Slice' - fg1d(gcoefs(:,3),xs));
figure(13); plot(image4Slice' - fg1d(gcoefs(:,4),xs));
figure(14); plot(image5Slice' - fg1d(gcoefs(:,5),xs));

figure(10); plot(image1Slice' - fg1d(gcoefs(:,1),xs),'color',[0 1 0]); hold on;
 plot(image2Slice' - fg1d(gcoefs(:,2),xs),'color',[0 0.8 0]);
 plot(image3Slice' - fg1d(gcoefs(:,3),xs),'color',[0 0.6 0]);
 plot(image4Slice' - fg1d(gcoefs(:,4),xs),'color',[0 0.4 0]);
 plot(image5Slice' - fg1d(gcoefs(:,5),xs),'color',[0 0.3 0]);

outputBinned1 = binMe(image1Slice' - fg1d(gcoefs(:,1),xs),xs,30);
outputBinned2 = binMe(image2Slice' - fg1d(gcoefs(:,2),xs),xs,30);
outputBinned3 = binMe(image3Slice' - fg1d(gcoefs(:,3),xs),xs,30);
outputBinned4 = binMe(image4Slice' - fg1d(gcoefs(:,4),xs),xs,30);
outputBinned5 = binMe(image5Slice' - fg1d(gcoefs(:,5),xs),xs,30);


figure(20); plot(outputBinned1(1,:),'color',[0 1 0]); hold on;
 plot(outputBinned2(1,:),'color',[0 0.8 0]);
 plot(outputBinned3(1,:),'color',[0 0.6 0]);
 plot(outputBinned4(1,:),'color',[0 0.4 0]);
 plot(outputBinned5(1,:),'color',[0 0.3 0]);

