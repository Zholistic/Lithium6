
TotalDatasets = 3;

directoryCell = cell(TotalDatasets,1);
datestring = cell(TotalDatasets,1);

directoryCell{1} = 'C:\Data\140318_833G_Crossover_Measurement_10us_Imagepulse_1_I\';
directoryCell{2} = 'C:\Data\140318_880G_Crossover_Measurement_10us_Imagepulse_1_Isat\';
directoryCell{3} = 'C:\Data\140317_997G_Crossover_Measurement_10us_Imagepulse_0_5_Isat\';

datestring{1} = '140318';
datestring{2} = '140318';
datestring{3} = '140317';

bins(1) = 18;
bins(2) = 18;
bins(3) = 18;


%directory = 'C:\Data\140409_2D_EOS_972G_10us_0_5isat_1us_10Isat_60Low_20high_int_CutFromOtherFolder\';
%directory = 'C:\Data\140318_880G_Crossover_Measurement_10us_Imagepulse_1_Isat\';
directory = 'C:\Data\140317_997G_Crossover_Measurement_10us_Imagepulse_0_5_Isat\';
date = '140317';
camera = 'top';
varstring = 'EvapVoltage';
pixelLength = 2.84e-6; %2.84 um topcam
massL6 = 9.988e-27; %9.988 x 10^27 kg
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
    %if(i>=length(fileLocList(:))-20)
    %    Isat = 135;
    %end
    
    %disp(num2str(Isat));
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imgArrayFresh(:,:,i) = atom2Image(:,:,1);
    
end

%Display every X image:
for i=1:length(imgArrayFresh(1,1,:))
    if(mod(i,15) == 0)       
        figure(i);
        imagesc(imgArrayFresh(:,:,i));        
    end
end

%Sort to re-order in terms of EvapVoltage:
ordImgArray = []; sortVarData = [];
[sortVarData, index] = sortrows(varData);
ordImgArray(:,:,:) = imgArrayFresh(:,:,index);

%Average the same EvapVoltage (X per EvapVoltage) without assuming continuity:
%Note: must be ordered...
addList = []; avgImgs = []; evapNums = [];
currEvap = sortVarData(1,1); k=1; j=0;
for i=1:length(ordImgArray(1,1,:))
    thisEvap = sortVarData(i,1);
    
    if thisEvap == currEvap
        addList(k) = i;
        k = k+1;
    else
        addList(k) = i;
        j=j+1;
        avgImgs(:,:,j) = centerAndAverage(ordImgArray(:,:,addList));
        evapNums(j) = currEvap;
        addList = []; k = 1;
        currEvap = thisEvap;
    end
end

%Now the fresh images are averaged over evap voltage. Center them and
%calculate ROI for each:
centeredImages = [];
centeredImages = centerImgArray(avgImgs);

%Cut out a better region:
%center of cloud at (118,86), region 40x40
ROIy = 31:141; %Top:Bot
ROIx = 63:173;

atomNums = []; cropImages = []; radProfiles = []; noiseValue = [];
cropImages = centeredImages(ROIy,ROIx,:);

%for i=1:length(cropImages(1,1,:))
%    atomNums(i) = sum(sum(cropImages(:,:,i)));
%end

%Now radially average the clouds...
for i=1:length(cropImages(1,1,:))
    radProfiles(:,:,i) = radAverageBigSquare(centeredImages(:,:,i));
    %radProfiles(:,:,i) = radAverage(cropImages(:,:,i));
    noiseValue(i) = mean(radProfiles(1,55:66,i));
    radProfiles(:,:,i) = (radProfiles(:,:,i) - noiseValue(i))*1.2; %1.2 is the 'real atom' scaling factor
end

for i=1:length(cropImages(1,1,:))
    atomNums(i) = sum(sum(cropImages(:,:,i) - noiseValue(i)));
end



coefsG = [];
%Fit a gaussian centered at zero:
for i=1:length(radProfiles(1,1,:))
    coefsG(:,i) = gausFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
end

%plot(radProfiles(2,:,30),radProfiles(1,:,30)); grid on; hold on; plot(radProfiles(2,:,30),fg(coefsG(:,30),radProfiles(2,:,30)));
%imgNo = 6; plot(radProfiles(2,:,imgNo),radProfiles(1,:,imgNo)); grid on; hold on; plot(radProfiles(2,:,imgNo),fg(coefsG(:,imgNo),radProfiles(2,:,imgNo)));

fg = @(p,x)(p(1).*exp((-1).*((x).^2) ./ (2.*p(2).^2))); 
for i=1:length(radProfiles(1,1,:))
    if(mod(i,5) == 0)       
        figure(i+1000);
        imgNo = i; 
        plot(radProfiles(2,:,imgNo),radProfiles(1,:,imgNo)); 
        grid on; hold on; 
        plot(radProfiles(2,:,imgNo),fg(coefsG(:,imgNo),radProfiles(2,:,imgNo)),'r');
        hold off;
    end
end

%Second Moment on the radial profiles:
SMom = []; ROIxSM = 1:60;
for i=1:length(radProfiles(1,1,:))
    SMom(i) = sum(radProfiles(1,ROIxSM,i).*(radProfiles(2,ROIxSM,i)).^2) / sum(radProfiles(1,ROIxSM,i));
end

%Region trick:
%plot(radProfiles(1,ROIxSM,50).*(radProfiles(2,ROIxSM,50)).^2)

sigmaX = sqrt(SMom);


%plot(atomNums,coefsG(2,:),'.');

disp('Computation post-denouement.');








