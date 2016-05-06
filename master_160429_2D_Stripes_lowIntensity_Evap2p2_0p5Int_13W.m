directory = 'C:\Data\160429_stripe_study_low_intensity_2D_Evap2p2_0p5Int_13Watt\';
%directory = 'C:\Data\151003_New_breath_test\';
%directory = 'C:\Data\150904_magnification_check3_top\';
%NOT ACTUALLY 724G. like 780 or 750, check lab book...
date = '160429';
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


%PCA Play%
if(0)
        imagesToReshape = [];
        imagesToReshape = imageArrayC(:,:,1:end); %single images
        %imagesToReshape = imagesAvgSet1;
        %imagesToReshape = imagesAvgSet2;
        %imagesToReshape = imagesAvgSet3;
        %imagesToReshape = imagesAvgSet4;
        
        numberOfImages = length(imagesToReshape(1,1,:));
        vectorLength = length(imagesToReshape(:,1,1))*length(imagesToReshape(1,:,1));
        
        %reshape images into vectors:
        vectorImages = []; vectorSum = zeros(vectorLength,1);
        for i=1:length(imagesToReshape(1,1,:))
            vectorImages(:,i) = reshape(imagesToReshape(:,:,i),[],1);
            vectorSum = vectorSum + vectorImages(:,i);
        end
        
        %Mean Image:
        meanImageVector = (1/numberOfImages).*vectorSum;
        undividedImage = reshape(vectorSum,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        meanImage = reshape(meanImageVector,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        
        normImages = [];
        for i=1:length(imagesToReshape(1,1,:))
            normImages(:,i) = vectorImages(:,i) - meanImageVector;
        end
        
        %N x p Bmatrix (N = number of images, p = pixels per image)
        Bmatrix = normImages;
        
        %[coeff, score, latent] = pca(ingredients)
        %each column of score corresponds to one principal component
        %latent stores the variances of the N principal components
        [Y01,P01,E01,tsquared,percentV] = pca(Bmatrix);
              
        
        %Reconstruct images using desired principle components
        eigenVectorBmatrix = []; imagesFromEVector = [];
        for i=1:numberOfImages
            eigenVectorBmatrix = P01(:,i);
            imagesFromEVector(:,:,i) = reshape(eigenVectorBmatrix,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        end
        
        if(0)
            for i=1:numberOfImages
                figure(i); imagesc(imagesFromEVector(:,:,i));
            end
        end
        
        
        %set of eigenmodes (principal components) large figure:
        figure(1);
        %set(gca,'XTickLabelMode', 'manual','XTickLabel', []);
        %set(gca,'YTickLabelMode', 'manual','YTickLabel', []);
        for i=1:numberOfImages
            subplot(floor(numberOfImages/5),5,i);
            imagesc(imagesFromEVector(:,:,i));
            set(gca,'XTickLabelMode', 'manual','XTickLabel', []);
            set(gca,'YTickLabelMode', 'manual','YTickLabel', []);         
        end

end     
%END PCA PLAY%


%minusTwenty = 1:10;

imageArray1 = imageArrayC(:,:,:);


imageArray1Avg = mean(imageArray1(:,:,:),3);
imageArray1Avg = centerAndAverage(imageArrayC(:,:,1:end));


fig80 = figure(80); imagesc(imageArray1Avg);



imageArray1AvgRot = imrotate(imageArray1Avg,-20);


fig80 = figure(8000); imagesc(imageArray1AvgRot);


xregionToAvg = 125:151;
yregionToAvg = 38:230;

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
    shiftByAvgs(i) = mean(radProfilesAvg(1,70:90,i));
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

