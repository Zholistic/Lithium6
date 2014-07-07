%Script for generating the equation of state compressibility vs pressure
%for the 60 low intensity\20 high intensity image sets. 

%directory = 'C:\Data\140616_2D_EOS_972G_No_Hold\';
directory = 'C:\Data\140704_2D_EOS_972G_15mv_evap\';
%directory = 'C:\Data\140411_2D_EOS_833G_10us_0_5isat_1us_10Isat_60Low_20high_int\';
%directory = 'C:\Data\140411_2D_EOS_880G_10us_0_5isat_1us_10Isat_60Low_20high_int\';
date = '140704';
camera = 'top';
varstring = 'Isat';
pixelLength = 2.84e-6; %2.84 um topcam
massL6 = 9.988e-27; %9.988 x 10^-27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
kB = 1.38e-23; %Boltzmanns constant m^2 kg s^-2 K^-1
imgArrayFresh = [];  lowIntRealAtomImg = [];
OD = 0; %optical density from SPE process function 1=OD, 0=WithSigma
close all;

%Get information from log file:
[fileLocList,varData] = generateFromLogfile(directory,date,varstring,camera);

%Build images from files:
for i=1:length(fileLocList(:))
    
    %Last 20 Images are at different Isat:
    if(varData(i,1) > 5) %12 Isat, high intensity images
        Isat = 135;
    end
    
    %disp(num2str(Isat));
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imgArrayFresh(:,:,i) = atom2Image(:,:,1);
    
end

%Display every X image:
for i=1:length(imgArrayFresh(1,1,:))
    if(mod(i,5) == 0)       
        figure(i);
        imagesc(imgArrayFresh(:,:,i));        
    end
end

lowIntImage = centerAndAverage(imgArrayFresh(:,:,51:end)); %works
highIntImage = centerAndAverage(imgArrayFresh(:,:,1:50));

lowIntImage = shiftToWingsZero(lowIntImage);
highIntImage = shiftToWingsZero(highIntImage);

%Center the images:
toCenter = []; centeredImages = [];
toCenter(:,:,2) = highIntImage;
toCenter(:,:,1) = lowIntImage;
centeredImages = centerImgArray(toCenter);

highIntImage = centeredImages(:,:,2);
lowIntImage = centeredImages(:,:,1);


figure(100)
imagesc(lowIntImage);
figure(101)
imagesc(highIntImage);


%scaledLowIntImage = lowIntImage.*toAtomScaleFactor;

%Specify region:
%ROIy = 50:110; ROIx = 70:145; %Good topcam region tightish
ROIy = 10:150; ROIx = 30:185; %Topcam looser region

highLowSpectrum = makeSpectrumHL(highIntImage(ROIy,ROIx), lowIntImage(ROIy,ROIx));
figure(102)
plot(highLowSpectrum(2,:),highLowSpectrum(1,:));
lowIntRealAtomImg = convertToRealAtomNumber(highIntImage(ROIy,ROIx), lowIntImage(ROIy,ROIx));
%lowIntRealAtomImg = convertToRealAtomNumber(highIntImage, lowIntImage);

NROISum = sum(sum(lowIntImage(ROIy,ROIx)));
%NROIAtoms = sum(sum(lowIntRealAtomImg(ROIy,ROIx)));
NROIAtoms = sum(sum(lowIntRealAtomImg));

%Center the images (again):
toCenter = []; centeredImages = [];
toCenter(:,:,2) = highIntImage(ROIy,ROIx);
toCenter(:,:,1) = lowIntImage(ROIy,ROIx);
toCenter(:,:,3) = lowIntRealAtomImg;
centeredImages = centerImgArray(toCenter);

highIntImage = centeredImages(:,:,2);
lowIntImage = centeredImages(:,:,1);
lowIntRealAtomImg = centeredImages(:,:,3);


%Profiles:
figure(206)
hold on;
plot(mean(lowIntImage(:,:),1));
plot(mean(lowIntRealAtomImg(:,:),1),'r');
plot(mean(highIntImage(:,:),1),'g');

figure(207)
hold on;
plot(mean(lowIntImage(:,:),2));
plot(mean(lowIntRealAtomImg(:,:),2),'r');
plot(mean(highIntImage(:,:),2),'g');


%Radial Profile:
radProfile = radAverage(lowIntRealAtomImg);
figure(203)
plot(radProfile(2,:),radProfile(1,:));