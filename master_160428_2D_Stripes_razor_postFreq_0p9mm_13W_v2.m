directory = 'C:\Data\160428_stripe_study_2D_aftertrapfreq_with_razorblade_at_0p9mm_13W_3v_v2\';
%directory = 'C:\Data\151003_New_breath_test\';
%directory = 'C:\Data\150904_magnification_check3_top\';
%NOT ACTUALLY 724G. like 780 or 750, check lab book...
date = '160428';
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
OD = 1; %optical density from SPE process function 1=OD, 0=WithSigma
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

%minusTwenty = 1:10;

imageArray1 = imageArrayC(:,:,:);


imageArray1Avg = mean(imageArray1(:,:,:),3);


fig80 = figure(80); imagesc(imageArray1Avg);


imageArray1AvgRot = imrotate(imageArray1Avg,-20);
%imageArray1AvgRot = imrotate(imageArray1Avg,+25);
%imageArray1AvgRot = imageArray1Avg;

fig80 = figure(80); imagesc(imageArray1AvgRot);

xregionToAvg = 125:151;
yregionToAvg = 38:230;
%xregionToAvg = 100:124;
%yregionToAvg = 5:180;

image1Slice = mean(imageArray1AvgRot(yregionToAvg,xregionToAvg),2);

%figure(1); subplot(5,1,1); plot(image1Slice);

gcoefsI = [];
gcoefsI(:,1) = gausFit1D(image1Slice);


%zeroing
image1Slice = image1Slice - mean(image1Slice(1:18));


gcoefs = [];
gcoefs(:,1) = gausFit1D(image1Slice);


fg1d = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));

xs = 1:length(image1Slice);

figure(1); plot(image1Slice); hold on; plot(fg1d(gcoefs(:,1),xs),'r');


figure(10); plot(image1Slice' - fg1d(gcoefs(:,1),xs));


figure(10); plot(image1Slice' - fg1d(gcoefs(:,1),xs),'color',[0 1 0]); hold on;


outputBinned1 = binMe(image1Slice' - fg1d(gcoefs(:,1),xs),xs,30);


figure(20); plot(outputBinned1(1,:),'color',[0 1 0]); hold on;


%--------------Temperature Fit: --------------%
if(0)

%Not centered before average:
imageArrayAvgs = imageArray1Avg;

%Centered before average:
imageArrayAvgs = centerAndAverage(imageArrayC(:,:,1:end));

%Radially averaged average profiles:
radProfilesAvg = []; radProfilesTAvg = []; centerAvgs = [];
disp('Radially averaging...');
for i=1:length(imageArrayAvgs(1,1,:))
    [radProfilesTAvg(:,:,i),centerAvgs(:,i)] = radAverageBigSquare(imageArrayAvgs(:,:,i));
    radProfilesAvg(:,:,i) = radProfilesTAvg(:,1:end-5,i);
end

%Shift wings of radial profiles to zero:
shiftByAvgs = [];
for i=1:length(radProfilesAvg(1,1,:))
    shiftByAvgs(i) = mean(radProfilesAvg(1,82:95,i))
    radProfilesAvg(1,:,i) = radProfilesAvg(1,:,i) - shiftByAvgs(i);
end

%Shift images by found factor:
imageArrayAvgsZerod = [];
for i=1:length(imageArrayAvgs(1,1,:))
    imageArrayAvgsZerod(:,:,i) = imageArrayAvgs(:,:,i) - shiftByAvgs(i);
end


if(0)
    for i=1:length(imageArrayAvgs(1,1,:))
        if(mod(i,1) == 0)
            figure(i);
            hold on; plot(radProfilesAvg(2,:,i),radProfilesAvg(1,:,i),'r'); line([0 140],[0 0],'color',[0 0 0]); hold off;
        end
    end
end



%make sure the y axis of rad profile is in 'pixels' not optical density. The 1.2 correction is
%done below
radproftoFityReal = radProfilesAvg(1,:,1)./(pixelLength^2);

%Fit Virial:
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
magfield = 832; %If you change this, also change a2d
omegaz = 5300 * 2 * pi; %Estimate at 5.3kHz, if you change this, also change a2d
omegar = 25 * 2 * pi; %Double check this weak trapping
az = sqrt(hbar/(massL6*omegaz));
pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
kpixelLength = (13e-6*(83/400));
radiusVector = 1:(length(radProfilesAvg(2,:,1)));
radiusVectorMeters = radiusVector.*kpixelLength;
potential = 0.5 .* massL6 .* omegar^2 .* radiusVectorMeters.^2;
potential_nk = (potential./kB).*(10^9);

%convert radProfile pixels to potential:
radProfilesPotential = 0.5 .* massL6 .* omegar^2 .*(radProfilesAvg(2,:,1).*pixelLength).^2;
%plot(radProfilesPotential,radProfilesAvg(1,:,2)./(pixelLength^2));
%radProfilesPotentialConv = radProfilesPotential.*10^30;
radProfilesDensityAdjusted = radproftoFityReal.*2.*1.2; %2 for spin states and 1.2 for correction factor

a0 = 5.29e-11; %o.0 A scaling factor...
a2d = 21420.8 * a0; %From mathematica spreadsheet TODO generate text file
eb = (hbar^2)/(a2d^2 * massL6);


%Scale x&y to H.O. dimensionless 
%Density SI = m^-2
radProfilesDensityConv = radProfilesDensityAdjusted./(massL6*omegaz / hbar);
%Potential SI = kg * s^-2 * m^2 // Energy units = kg * m^2 * s^-2
radProfilesPotentialConv = radProfilesPotential./((hbar*omegaz)/2); 

figure(222); plot(radProfilesPotentialConv(1:70),radProfilesDensityConv(1:70),'.');

fitrange = 35:50; %as soon as it goes below zero it's going to bug out, 10 to 15 points is sufficient

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
    
    figure(3); plot(virialFitFuncs{13}); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.');
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

end
