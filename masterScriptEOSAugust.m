%Script for generating the equation of state compressibility vs pressure
%for the 60 low intensity\20 high intensity image sets. 

%directory = 'C:\Data\140814_2D_EOS_972G_2D_evap_1p25v_10us_0p5isat_1us_10isat\';
%directory = 'C:\Data\140814_2D_EOS_972G_2D_evap_1p25v_10us_0p5isat_1us_10isat_4sec_hold\';
directory = 'C:\Data\140819_2D_EOS_880G_2D_evap_1p2v_10us_0p3isat_1us_10isat\';
date = '140819';
camera = 'top';
varstring = 'Isat';
%varstring = 'isat'; %972G Datasets
pixelLength = 2.84e-6; %2.84 um topcam
massL6 = 1e-26; %9.988 x 10^-27 kg
hbar = 1.055e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
kB = 1.38e-23; %Boltzmanns constant m^2 kg s^-2 K^-1
imgArrayFresh = [];  lowIntRealAtomImg = [];
OD = 0; %optical density from SPE process function 1=OD, 0=WithSigma
close all; varData = [];

%Get information from log file:
[fileLocList,varData] = generateFromLogfile(directory,date,varstring,camera);

%Build images from files:
for i=1:length(fileLocList(:))
    Isat = 10*135;
    %Last 20 Images are at different Isat:
    if(varData(i,1) > 0.6)
        Isat = 135;
    end
    
    %disp(num2str(Isat));
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imgArrayFresh(:,:,i) = atom2Image(:,:,1);
    
end

%Display every X image:
if(0)
for i=1:length(imgArrayFresh(1,1,:))
    if(mod(i,5) == 0)       
        figure(i);
        imagesc(imgArrayFresh(:,:,i));        
    end
end
end

lowIntImage = centerAndAverage(imgArrayFresh(:,:,1:200));
highIntImage = centerAndAverage(imgArrayFresh(:,:,201:end));

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

spectrumFunc1 = makeSpectrumHL(highIntImage(ROIy,ROIx),lowIntImage(ROIy,ROIx));

highLowSpectrum = makeSpectrumHL(highIntImage(ROIy,ROIx), lowIntImage(ROIy,ROIx));
if(0)
    figure(102)
    plot(highLowSpectrum(2,:),highLowSpectrum(1,:));
end

scaleFactor = mean(spectrumFunc1(1,40:80));

lowIntRealAtomImg = lowIntImage(ROIy,ROIx).*scaleFactor;

%lowIntRealAtomImg = convertToRealAtomNumber(highIntImage(ROIy,ROIx), lowIntImage(ROIy,ROIx));
%lowIntRealAtomImg = convertToRealAtomNumber(highIntImage, lowIntImage);

NROISum = sum(sum(lowIntImage(ROIy,ROIx)));
%NROIAtoms = sum(sum(lowIntRealAtomImg(ROIy,ROIx)));
NROIAtoms = sum(sum(lowIntRealAtomImg));

if(0) %%%%%%----------%----------%
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

%just plot a strip in the middle, not an average:
figure(202)
hold on;
plot(lowIntImage(70,:));
plot(lowIntRealAtomImg(70,:),'r');
plot(highIntImage(70,:),'g');
end %%%%%%----------%----------%


%Radial Profile:
radProfile = []; center = [];
%radProfile = radAverage(lowIntRealAtomImg);
[radProfile,center] = radAverageBigSquare(lowIntRealAtomImg);
figure(203)
plot(radProfile(2,:),radProfile(1,:)); %Atoms per pixel...

%Shift wings to zero:
offset = mean(radProfile(1,60:80));
radProfile(1,:) = radProfile(1,:) - offset;

%Polylog fit to profile:
gcoefsPolyLog1 = [];
fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*x.^2)./(p(3).^2))));
gcoefsPolyLog1 = polyLog1FitHalf1D(radProfile(1,:),radProfile(2,:));

if(0)
    figure(203);
    plot(radProfile(2,:),radProfile(1,:));
    hold on;
    plot(fgp(gcoefsPolyLog1,1:100),'r');
end

radProfilePolylog(1,:) = fgp(gcoefsPolyLog1,1:length(radProfile(1,:)));
radProfilePolylog(2,:) = radProfile(2,:);

%Physical calculations:
TonTF = 1/(log(1+exp(gcoefsPolyLog1(2)/gcoefsPolyLog1(3)^2)));
Temp = (2*pi*((hbar)^2)*(gcoefsPolyLog1(1)/(pixelLength^2)))/(massL6*kB);
%pixelLengthAvg = ((2^(1/2))*pixelLength + pixelLength)/2; %Radial average pixel length!
omegaR = 2*pi*23.7; %Radial trapping frequency (23.7 Hz)
radLength = radProfile(2,:).*pixelLength; %Pixel length in the hypotenuse?
radPotential = (massL6)*(0.5)*(((omegaR).*(radLength)).^2);
radTemp = (radPotential./kB)./(10^-9);
%radProfile(1,:) = radProfile(1,:)./(10^-12); %convert from um^-2 to m^-2 area
%convert from pixel number to density at that pixel
%radProfile(1,:) = radProfile(1,:)./(10^-12); %convert from pixel area um^-2 to m^-2 area
radProfile(1,:) = radProfile(1,:)./(pixelLength^2);

%Polylog 'ideal' EOS physical data:
radLengthPolylog = radProfilePolylog(2,:).*pixelLength;
radPotentialPolylog = (massL6)*(0.5)*(((omegaR).*(radLengthPolylog)).^2);
radTempPolylog = (radPotentialPolylog./kB)./(10^-9);
radProfilePolylog(1,:) = radProfilePolylog(1,:)./(pixelLength^2);

%Ideal calculations:
PIdeal = []; KIdeal = [];
%Center pixel atom number with N/A density...
%radProfile is in atoms/m^2
K0 = (massL6/(pi*(hbar^2))).*(1./((radProfile(1,1)*(10^-12))/(pixelLength^2).^2)); 
%P0 = (pi*(hbar^2)/(2*massL6)).*(((radProfile(1,1)*(10^-12))/(pixelLength^2)).^2);
P0 = (pi*(hbar^2)/(2*massL6))*((radProfile(1,1))^2);
%P0D = (pi*(hbar^2)/(2*massL6))*((radDensity(1,1))^2);
P = trapz(radPotential(5:end-20),radProfile(1,5:end-20));
%P = sum(radPotential(5:end-20))*trapz(radProfile(1,5:end-20));
%PD = sum(radDensity(2,1:end-20))*trapz(radDensity(1,1:end-20));

if(0)
figure(204)
plot(radTemp,radProfile(1,:),'.'); %density vs V(r) in nK

figure(205)
plot(radPotential(1:end-20),radProfile(1,1:end-20),'.'); %density vs V(r) in nK
end

%Re-bin the profile vs potential:
nbins = 100; binMean = []; radialDensity = []; binStd = []; binEdges=[];
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
radialDensity(2,i) = min(radPotential) - binLengths(1) +sum(binLengths(1:i));
end

P0D = (pi*(hbar^2)/(2*massL6))*((radialDensity(1,1))^2);
PD = trapz(radialDensity(2,1:end-20),radialDensity(1,1:end-20));

figure(750)
plot(radialDensity(2,:),radialDensity(1,:),'.');


%X_1:
X1 = []; X1BD = []; X1Polylog = []; %BD = Binned Density
%X1 = (-1)*(hbar^2 / massL6).*diff(radProfile(1,:));
%X1 = (-1)*diff(radProfile(1,:));
%X1 = (-1)*(hbar^2 / massL6).*gradient(radProfile(1,:),radPotential);
X1 = (-1)*((hbar^2) / massL6).*gradient(radProfile(1,8:end-20),radPotential(8:end-20));
X1Polylog = (-1)*((hbar^2) / massL6).*gradient(radProfilePolylog(1,1:end),radPotentialPolylog(1:end));
%diffF = diff(radialDensity(1,1:end-20))./diff(radialDensity(2,1:end-20));
X1BD = (-1)*((hbar^2) / massL6).*gradient(radialDensity(1,1:end-22),radialDensity(2,1:end-22)); %gradient(y,x)

%X_(-1):
Xm1 = []; Xm1BD = []; %BD = Binned Density
Xm1Polylog = [];
%for i=1:length(radProfile(1,:))
for i=8:(length(radProfile(1,1:end-20))-1)
  
    Xm1(i-7) = (massL6 / (hbar^2))*(1/(radProfile(1,i)^2)).*trapz(radPotential(i:end-20),radProfile(1,i:end-20)); %may be backwards
    %Xm1(i) = (1/(radProfile(1,i)^2)).*trapz(radProfile(1,i:end));
    
end

for i=1:(length(radialDensity(1,1:end-22)))
  
    Xm1BD(i) = (massL6 / (hbar^2))*(1/(radialDensity(1,i)^2))*trapz(radialDensity(2,i:end-22),radialDensity(1,i:end-22),2); %trapz(x)*spacing
    
end

for i=1:(length(radProfilePolylog(1,1:end)))
  
    Xm1Polylog(i) = (massL6 / (hbar^2))*(1/(radProfilePolylog(1,i)^2))*trapz(radPotentialPolylog(i:end),radProfilePolylog(1,i:end),2); %trapz(x)*spacing
    
end

%Pressure and Compressibility:
PonP0 = (2/pi).*Xm1;
KonK0L = pi.*X1;
KonK0 = KonK0L(1:end-1);

%For Binned Density:
PonP0BD = (2/pi).*Xm1BD(1:end-2);
KonK0LBD = pi.*X1BD(1:end-2);
KonK0BD = KonK0LBD(1:end);

%For Polylog Function:
PonP0PL = (2/pi).*Xm1Polylog;
KonK0LPL = pi.*X1Polylog;
KonK0PL = KonK0LPL(1:end-1);

if(0)
 figure(250);
 plot(Xm1BD);
 figure(251);
 plot(X1BD); 
 plot(PonP0PL(1:48),KonK0LPL(1:48),'.')
end



%Bin both the comp and pressure:
%KonK0Binned = []; PonP0Binned = [];
%KonK0Binned = binMe(KonK0BD(1:48),1:length(KonK0BD(1:48)),50);
%PonP0Binned = binMe(PonP0BD(1:48),1:length(PonP0BD(1:48)),50);


%Output plots:
%figure(500)
%plot(radTemp,radProfile(1,:).*(10^-12));
%title('Density vs Potential at 880G'); xlabel('Potential V(r) in nK'); ylabel('Density n in um^2');
    
%figure(501)
%plot(radLength,radProfile(1,:).*(10^-12)); grid on;
%title('Density vs Width at 880G'); xlabel('Width in m'); ylabel('Density n in um^2');

figure(502)
plot(KonK0,PonP0,'.');

figure(512)
%plot(KonK0BD,PonP0BD,'.');
plot(PonP0BD,KonK0BD,'.');
title('Binned on Density Comp vs Pressure');

%figure(513)
%plot(KonK0Binned(1,:),PonP0Binned(1,:),'.');
%title('Binned on each of comp and pressure Comp vs Pressure');

%Individual Bin for the comp/pressure bins:
%KonK0IBb = []; PonP0IBb = []; KvsPb = [];
%onK0IBb = binMe(KonK0Binned(1,:),1:length(KonK0Binned(1,:)),60);
%PonP0IBb = binMe(PonP0Binned(1,:),1:length(PonP0Binned(1,:)),60);

%loop to cut out low P/P0 values (less than 'cutoffPP0'):
%cutoffPP0 = 0.8;
%for i=lengthPonP0BD(1,:)
%end


KvsPb = [];
%KvsPb = binMe(PonP0BD(1,:),KonK0BD(1,:),25);
KvsPb = binMe(KonK0BD,PonP0BD,50);

figure(504)
errorbar(KvsPb(2,:),KvsPb(1,:),KvsPb(3,:)./2,'.');
%plot(KvsPb(2,:),KvsPb(1,:),'.');
%errorbar(KvsPb(1,:),KvsPb(2,:),KvsPb(3,:)./2,'.');


%%%%Plot with line overlay:
if(0)
figure(504); hold on;
%errorbar(KvsP833(2,:),KvsP833(1,:),KvsP833(3,:)./2,'.');
plot(KvsP833(2,:),KvsP833(1,:));
errorbarxy(KvsP833(2,:),KvsP833(1,:),zeros(length(KvsP833(1,:))),KvsP833(3,:)./2,zeros(length(KvsP833(1,:))),KvsP833(3,:)./2,'.');
%errorbar(KvsP880(2,:),KvsP880(1,:),KvsP880(3,:)./2,'.r');
hold on; plot(KvsP880(2,:),KvsP880(1,:),'r');
errorbarxy(KvsP880(2,:),KvsP880(1,:),zeros(length(KvsP880(1,:))),KvsP880(3,:)./2,zeros(length(KvsP880(1,:))),KvsP880(3,:)./2,'.');
%errorbar(KvsP972(2,:),KvsP972(1,:),KvsP972(3,:)./2,'.g');
hold on; plot(KvsP972(2,:),KvsP972(1,:),'g');
errorbarxy(KvsP972(2,:),KvsP972(1,:),zeros(length(KvsP972(1,:))),KvsP972(3,:)./2,zeros(length(KvsP972(1,:))),KvsP972(3,:)./2,'.');
end
%%%%

if(0)
figure(504); hold on;
errorbarxy(KvsP833(2,:),KvsP833(1,:),zeros(length(KvsP833(1,:))),KvsP833(3,:)./2,zeros(length(KvsP833(1,:))),KvsP833(3,:)./2,'.');
errorbarxy(KvsP880(2,:),KvsP880(1,:),zeros(length(KvsP880(1,:))),KvsP880(3,:)./2,zeros(length(KvsP880(1,:))),KvsP880(3,:)./2,'.');
errorbarxy(KvsP972(2,:),KvsP972(1,:),zeros(length(KvsP972(1,:))),KvsP972(3,:)./2,zeros(length(KvsP972(1,:))),KvsP972(3,:)./2,'.');
end

%%plot(KonK0(1:45),PonP0(1:45),'.');

% figure(1000); plot(KonK0_972,PonP0_972,'.'); hold on; plot(KonK0_880,PonP0_880,'^'); plot(KonK0_833,PonP0_833,'s');

%PonP0 = PonP0Binned;
%KonK0 = KonK0Binned;

%Individual Bin:
%KonK0IB = []; PonP0IB = []; KvsP = [];
%KonK0IB = binMe(KonK0,1:length(KonK0),60);
%PonP0IB = binMe(PonP0,1:length(PonP0),60);

%KvsP = binMe(KonK0IB(1,:),PonP0IB(1,:),25);

%figure(503)
%plot(KvsP(1,:),KvsP(2,:),'.');

% figure(1000); plot(KvsP_972(2,:),KvsP_972(1,:),'.'); hold on; plot(KvsP_880(2,:),KvsP_880(1,:),'^'); plot(KvsP_833(2,:),KvsP_833(1,:),'s');
%plot(Xm1BD,X1BD,'.'); axis([0 12 0 4]);










