%Script for generating the equation of state compressibility vs pressure
%for the 60 low intensity\20 high intensity image sets. 

%directory = 'C:\Data\140616_2D_EOS_972G_No_Hold\';
directory = 'C:\Data\140704_2D_EOS_972G_15mv_evap\';
%directory = 'C:\Data\140707_2D_EOS_972G_15mv_evap_2sechold\';
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
    
    Isat = 10*135;
    
    %High Intensity images different isat:
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

ROIy = 10:150; ROIx = 30:185; %Topcam looser region
rawNums = [];
for i=1:length(imgArrayFresh(1,1,:))
    rawNums(i) = sum(sum(imgArrayFresh(ROIy,ROIx,i)));
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
radProfile = radAverageBigSquare(lowIntRealAtomImg);
figure(203)
plot(radProfile(2,:),radProfile(1,:));

%Physical calculations:
%pixelLengthAvg = ((2^(1/2))*pixelLength + pixelLength)/2; %Radial average pixel length!
omegaR = 2*pi*23.7; %Radial trapping frequency (23.7 Hz)
radLength = radProfile(2,:).*pixelLength; %Pixel length in the hypotenuse?
radPotential = (massL6)*(0.5)*(((omegaR).*(radLength)).^2);
radTemp = (radPotential./kB)./(10^-9);
%radProfile(1,:) = radProfile(1,:)./(10^-12); %convert from um^-2 to m^-2 area
%convert from pixel number to density at that pixel
radProfile(1,:) = radProfile(1,:)./(10^-12); %convert from um^-2 to m^-2 area


%Ideal calculations:
PIdeal = []; KIdeal = [];
%Center pixel atom number with N/A density...
%radProfile is in atoms/m^2
K0 = (massL6/(pi*(hbar^2))).*(1./((radProfile(1,1)*(10^-12))/(pixelLength^2).^2)); 
%P0 = (pi*(hbar^2)/(2*massL6)).*(((radProfile(1,1)*(10^-12))/(pixelLength^2)).^2);
P0 = (pi*(hbar^2)/(2*massL6))*((radProfile(1,1))^2);
%P0D = (pi*(hbar^2)/(2*massL6))*((radDensity(1,1))^2);
%P = trapz(radPotential(5:end-20),radProfile(1,5:end-20));
P = sum(radPotential(5:end-20))*trapz(radProfile(1,5:end-20));
%PD = sum(radDensity(2,1:end-20))*trapz(radDensity(1,1:end-20));


figure(204)
plot(radTemp,radProfile(1,:),'.'); %density vs V(r) in nK

figure(205)
plot(radPotential(1:end-20),radProfile(1,1:end-20),'.'); %density vs V(r) in nK

%Re-bin the profile vs potential:
nbins = 160; binMean = []; radialDensity = []; binStd = []; binEdges=[];
binEdges = linspace(min(radPotential),max(radPotential),nbins+1);
binLength = (max(radPotential)-min(radPotential))/(nbins+1);
[h,whichBin] = histc(radPotential, binEdges);

j=1; binCount = 1;
for i = 1:nbins
    flagBinMembers = (whichBin == i);
    binMembers     = radProfile(1,flagBinMembers);
    %what bins to exclude? NaN and really small values (less than 1)
    if ~isnan(mean(binMembers)) && (mean(binMembers) > 1)
        binMean(j)     = mean(binMembers);
        binStd(j)      = std(binMembers);
        binLengths(j) = binCount*binLength;
        j=j+1;
        binCount = 1;
    else
        binCount = binCount +1;
    end
end


for i=1:length(binMean)
radialDensity(1,i) = binMean(i);
radialDensity(2,i) = sum(binLengths(1:i));
end

P0D = (pi*(hbar^2)/(2*massL6))*((radialDensity(1,1))^2);
PD = trapz(radialDensity(2,1:end-20),radialDensity(1,1:end-20));

figure(750)
plot(radialDensity(2,:),radialDensity(1,:),'.');

%X_1:
X1 = []; X1BD = []; %BD = Binned Density
X1 = (-1)*((hbar^2) / massL6).*gradient(radProfile(1,8:end-20),radPotential(8:end-20));
X1BD = (-1)*((hbar^2) / massL6).*gradient(radialDensity(1,1:end-25),radialDensity(2,1:end-25)); %gradient(y,x)

%X_(-1):
Xm1 = []; Xm1BD = []; %BD = Binned Density
for i=8:(length(radProfile(1,1:end-20))-1)
  
    Xm1(i-7) = (massL6 / (hbar^2))*(1/(radProfile(1,i)^2)).*trapz(radPotential(i:end-20),radProfile(1,i:end-20)); %may be backwards
    %Xm1(i) = (1/(radProfile(1,i)^2)).*trapz(radProfile(1,i:end));
    
end

for i=1:(length(radialDensity(1,1:end-25)))
  
    Xm1BD(i) = (massL6 / (hbar^2))*(1/(radialDensity(1,i)^2))*trapz(radialDensity(2,i:end-25),radialDensity(1,i:end-25),2); %trapz(x)*spacing
    
end

%Pressure and Compressibility:
PonP0 = (2/pi).*Xm1;
KonK0L = pi.*X1;
KonK0 = KonK0L(1:end-1);

%For Binned Density:
PonP0BD = (2/pi).*Xm1BD;
KonK0LBD = pi.*X1BD;
KonK0BD = KonK0LBD(1:end);



figure(512)
plot(KonK0BD,PonP0BD,'.');
title('Binned on Density Comp vs Pressure');


KvsPb = [];
KvsPb = binMe(PonP0BD(1,:),KonK0BD(1,:),12);

figure(504)
%plot(KvsPb(2,:),KvsPb(1,:),'.');
errorbar(KvsPb(2,:),KvsPb(1,:),KvsPb(3,:),'.');





