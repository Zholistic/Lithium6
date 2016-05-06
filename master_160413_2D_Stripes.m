directory = 'C:\Data\160413_2D_10W_Stripes\';
%directory = 'C:\Data\151003_New_breath_test\';
%directory = 'C:\Data\150904_magnification_check3_top\';
%NOT ACTUALLY 724G. like 780 or 750, check lab book...
date = '160413';
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

origpos = 9:18;
cwsome = 19:28;
cwfull = 29:38;
acwfull = 39:48;

imageArray1 = imageArrayC(:,:,origpos);
imageArray2 = imageArrayC(:,:,cwsome);
imageArray3 = imageArrayC(:,:,cwfull);
imageArray4 = imageArrayC(:,:,acwfull);

imageArray1Avg = mean(imageArray1(:,:,:),3);
imageArray2Avg = mean(imageArray2(:,:,:),3);
imageArray3Avg = mean(imageArray3(:,:,:),3);
imageArray4Avg = mean(imageArray4(:,:,:),3);

figure(1); imagesc(imageArray1Avg);
figure(2); imagesc(imageArray2Avg);
figure(3); imagesc(imageArray3Avg);
figure(4); imagesc(imageArray4Avg);
