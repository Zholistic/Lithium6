%Script for generating the equation of state compressibility vs pressure
%for the 60 low intensity\20 high intensity image sets. 

directory = 'C:\Data\140409_2D_EOS_972G_10us_0_5isat_1us_10Isat_60Low_20high_int_CutFromOtherFolder\';
%directory = 'C:\Data\140411_2D_EOS_833G_10us_0_5isat_1us_10Isat_60Low_20high_int\';
%directory = 'C:\Data\140411_2D_EOS_880G_10us_0_5isat_1us_10Isat_60Low_20high_int\';
date = '140409';
camera = 'top';
varstring = 'BField';
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
    if(i>=length(fileLocList(:))-20)
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

lowIntImage = centerAndAverage(imgArrayFresh(:,:,1:60));
highIntImage = centerAndAverage(imgArrayFresh(:,:,61:end));

%the centered low intensity images:
lowIntImages = []; lowIntRealAtomImages = [];
lowIntImages = centerImgArray(imgArrayFresh(:,:,1:60)); 


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

%Operations on the image set:
NROIAtomsEachImage = []; radProfilesEI = []; NROIAtomsEachImageBeforeConvert = [];
for i=1:length(lowIntImages(1,1,:))
    lowIntImages(:,:,i) = shiftToWingsZero(lowIntImages(:,:,i));
    lowIntRealAtomImages(:,:,i) = convertToRealAtomNumber(highIntImage(ROIy,ROIx), lowIntImages(ROIy,ROIx,i));
    NROIAtomsEachImage(i) = sum(sum(lowIntRealAtomImages(:,:,i)));
    NROIAtomsEachImageBeforeConvert(i) = sum(sum(lowIntImages(ROIy,ROIx,i)));
    radProfilesEI(:,:,i) = radAverage(lowIntRealAtomImages(:,:,i));
end
    


NROISum = sum(sum(lowIntImage(ROIy,ROIx)));
%NROIAtoms = sum(sum(lowIntRealAtomImg(ROIy,ROIx)));
NROIAtoms = sum(sum(lowIntRealAtomImg));

%Center the images (again):
toCenter = []; centeredImages = [];
toCenter(:,:,2) = highIntImage(ROIy,ROIx);
toCenter(:,:,1) = lowIntImage(ROIy,ROIx);
toCenter(:,:,3) = lowIntRealAtomImg;
centeredImages = centerImgArray(toCenter);

highIntImage = [];
highIntImage = centeredImages(:,:,2);
lowIntImage = centeredImages(:,:,1);
lowIntRealAtomImg = centeredImages(:,:,3);


%Profiles:
%200 is y, 201 is x:
figure(200)
hold on;
plot(mean(lowIntImage(:,:),1));
plot(mean(lowIntRealAtomImg(:,:),1),'r');
plot(mean(highIntImage(:,:),1),'g');

figure(201)
hold on;
plot(mean(lowIntImage(:,:),2));
plot(mean(lowIntRealAtomImg(:,:),2),'r');
plot(mean(highIntImage(:,:),2),'g');


%Radial Profile:
%Radially average each image, not the set:
for i=length(lowIntRealAtomImages(1,1,:))
end


radProfile = radAverage(lowIntRealAtomImg);
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
K0 = (massL6/(pi*(hbar^2))).*(1./((radProfile(1,1)*(10^-12))/(pixelLength^2).^2)); 
P0 = (pi*(hbar^2)/(2*massL6)).*(((radProfile(1,1)*(10^-12))/(pixelLength^2)).^2);

figure(204)
plot(radTemp,radProfile(1,:),'.'); %density vs V(r) in nK

figure(205)
plot(radPotential,radProfile(1,:),'.');


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

figure(750)
plot(radialDensity(2,:),radialDensity(1,:),'.');


%X_1:
X1 = []; X1BD = []; %BD = Binned Density
%X1 = (-1)*(hbar^2 / massL6).*diff(radProfile(1,:));
%X1 = (-1)*diff(radProfile(1,:));
%X1 = (-1)*(hbar^2 / massL6).*gradient(radProfile(1,:),radPotential);
X1 = (-1)*((hbar^2) / massL6).*gradient(radProfile(1,10:end-30),radPotential(10:end-30));
X1BD = (-1)*((hbar^2) / massL6).*gradient(radialDensity(1,1:end-25),radialDensity(2,1:end-25)); %gradient(y,x)

%X_(-1):
Xm1 = []; Xm1BD = []; %BD = Binned Density
%for i=1:length(radProfile(1,:))
for i=10:(length(radProfile(1,1:end-30))-1)
  
    Xm1(i-9) = (massL6 / (hbar^2))*(1/((radProfile(1,i))^2)).*trapz(radPotential(i:end-30),radProfile(1,i:end-30)); %may be backwards
    %Xm1(i) = (1/(radProfile(1,i)^2)).*trapz(radProfile(1,i:end));
    
end

for i=1:(length(radialDensity(1,1:end-25)))
  
    Xm1BD(i) = (massL6 / (hbar^2))*(1/(radialDensity(1,i)^2)).*trapz(radialDensity(2,i:end-25),radialDensity(1,i:end-25),2); %trapz(x,y)
    
end


%Pressure and Compressibility:
PonP0 = (2/pi).*Xm1;
KonK0L = pi.*X1;
KonK0 = KonK0L(1:end-1);

%For Binned Density:
PonP0BD = (2/pi).*Xm1BD;
KonK0LBD = pi.*X1BD;
KonK0BD = KonK0LBD(1:end);

%Bin both the comp and pressure:
KonK0Binned = []; PonP0Binned = [];
KonK0Binned = binMe(KonK0BD(1:48),1:length(KonK0BD(1:48)),50);
PonP0Binned = binMe(PonP0BD(1:48),1:length(PonP0BD(1:48)),50);


%Output plots:
figure(500)
plot(radTemp,radProfile(1,:).*(10^-12));
title('Density vs Potential at 880G'); xlabel('Potential V(r) in nK'); ylabel('Density n in um^(-2)');
    
figure(501)
plot(radLength,radProfile(1,:).*(10^-12)); grid on;
title('Density vs Width at 880G'); xlabel('Width in m'); ylabel('Density n in um^(-2)');

figure(502)
plot(PonP0,KonK0,'.');

figure(512)
plot(PonP0BD,KonK0BD,'.');
title('Binned on Density Comp vs Pressure');

figure(513)
plot(PonP0Binned(1,:),KonK0Binned(1,:),'.');
title('Binned on each of comp and pressure Comp vs Pressure');

%%plot(KonK0(1:45),PonP0(1:45),'.');

% figure(1000); plot(KonK0_972,PonP0_972,'.'); hold on; plot(KonK0_880,PonP0_880,'^'); plot(KonK0_833,PonP0_833,'s');

%PonP0 = PonP0Binned;
%KonK0 = KonK0Binned;

%Individual Bin:
KonK0IB = []; PonP0IB = [];
KonK0IB = binMe(KonK0,1:length(KonK0),60);
PonP0IB = binMe(PonP0,1:length(PonP0),60);

KvsP = binMe(KonK0IB(1,:),PonP0IB(1,:),100);

figure(503)
plot(KvsP(2,:),KvsP(1,:),'.');

%figure(504)
%plot(KvsP(1,:),KvsP(2,:),'.');

% figure(1000); plot(KvsP_972(2,:),KvsP_972(1,:),'.'); hold on; plot(KvsP_880(2,:),KvsP_880(1,:),'^'); plot(KvsP_833(2,:),KvsP_833(1,:),'s');
%plot(Xm1BD,X1BD,'.'); axis([0 12 0 4]);










