%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TOP
directory = 'C:\Data\140911_atomnumber_vs_magneticevapvoltage_topcam\';
date = '140911';
camera = 'top';
varstring = 'motfet';
%varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %2.84 um topcam, 
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
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

varDataTop = [];

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
[fileLocList,varDataTop] = generateFromLogfile(directory,date,varstring,camera);


imageArrayTop = [];
for i=1:length(fileLocList(:))
    
    Isat = 10*135;
    
    %High Intensity images different isat:
    if(0) %NO High intensity data for this run
        Isat = 135;  
    end
    
    %disp(num2str(Isat));
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArrayTop(:,:,i) = atom2Image(:,:,1);
    
end

imageArrayC = []; imageArrayTC = []; imageArrayHighIntensityC = [];
imageArrayHighIntensityTC = [];
ROIx = 1:length(imageArrayTop(1,:,1));
ROIy = 1:length(imageArrayTop(:,1,1));
CrossROIy = 10:150; 
CrossROIx = 30:185; %The cross is inside the region specified above.
TightROIx = 30:185;
TightROIy = 10:150;
%Split into high and low intensity arrays
imageArrayC = imageArrayTop(ROIy,ROIx,:);
%imageArrayHighIntensityC = imageArray(ROIy,ROIx,381:end); 
imageArrayTC = imageArrayTop(TightROIy,TightROIx,:);
%imageArrayHighIntensityTC = imageArray(TightROIy,TightROIx,381:end); 

finalScaleFactor = 1.3;

for i=1:length(imageArrayC(1,1,:))
    imageArrayC(:,:,i) = imageArrayC(:,:,i).*finalScaleFactor;
end

%%%%%Atom numbers:
%Tight ROI array:
pixelCountsTop = [];
for i=1:length(imageArrayC(1,1,:))
    pixelCountsTop(i) = sum(sum(imageArrayC(:,:,i)));
end

%%%%%%%%%%%%%%%%SIDE
directory = 'C:\Data\140911_atomnumber_vs_magneticevapvoltage_sidecam\';
date = '140911';
camera = 'sidecam';
varstring = 'motfet';

varDataSide = [];

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
[fileLocList,varDataSide] = generateFromLogfile(directory,date,varstring,camera);


imageArraySide = [];
%Pull images:
for i=1:length(fileLocList)
    imageArraySide(:,:,i) = PullFTS(fileLocList{i},raw);
end

%Side pixel counts are from the logfile!
pixelCountsSide = varDataSide(:,5);

for i=1:length(pixelCountsSide)
ratioST(i) = pixelCountsSide(i)/pixelCountsTop(i);
%ratioSTLog(i) = varDataSide(i,5)/varDataTop(i,5);
end

ratioTS = 1./ratioST;

figure(1);
plot(varDataSide(:,1),varDataSide(:,5));
hold on;
plot(varDataTop(:,1).*100,varDataTop(:,5),'r');
plot(varDataTop(:,1).*100,pixelCountsTop,'g');
hold off;

figure(2);
plot(varDataSide(:,1),ratioTS);

finalCorrection = mean(ratioTS)




