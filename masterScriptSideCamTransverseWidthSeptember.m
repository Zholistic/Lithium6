%directory = 'C:\Data\140904_transversewidth_600usTOF_Isat1e8_alpha0.6_atomnumber_5000\';
%directory = 'C:\Data\140904_transversewidth_600usTOF_Isat1e8_alpha0.6_atomnumber_10000\';
directory = 'C:\Data\140907_2D_transversewidth_10k_atoms_750ms_ramp\';
%directory = 'C:\Data\140909_transversewidth_5k_atoms_750ms_ramp\';
%directory = 'C:\Data\140915_transversewidth_13p5k_atoms_750ms_ramp\';
date = '140904';
date = '140907';
%date = '140909';
%date = '140915';
camera = 'sidecam';
varstring = 'mag_field'; %140907 (and previous?)
%varstring = 'magnetic_field';
atomnum = 10000;
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
[fileLocList,varData] = generateFromLogfile(directory,date,varstring,camera);

imageArray = [];
%Pull images:
for i=1:length(fileLocList)
    imageArray(:,:,i) = PullFTS(fileLocList{i},raw);
end


%Crop images:
imageArrayC = []; imageArrayTC = [];
ROIx = 450:730;
ROIy = 500:675;
%The cross is inside the region specified above. 
%Can improve the error by tightening the Cross ROI up
%for that particular dataset. 
CrossROIx = 75:200; 
CrossROIy = 65:100;

if(strcmp(date,'140907'))
    disp('140907 Data');
    %CrossROIx = 35:255; %The cross is inside the region specified above.
    %CrossROIy = 40:130;
end
    
TightROIx = 510:670;
TightROIy = 575:595;
imageArrayC = imageArray(ROIy,ROIx,:);
imageArrayTC = imageArray(TightROIy,TightROIx,:);


%Display every X image:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,5) == 0)       
        figure(i);
        imagesc(imageArrayC(:,:,i));        
    end
end
end


 
%%%%%Fit Gaussians:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
for i=1:length(imageArrayC(1,1,:))   
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    
    %Shift wings to zero:
    shiftFactor = (gcoefsXi(4,i)+gcoefsYi(4,i))/2;
    imageArrayC(:,:,i) = imageArrayC(:,:,i) - shiftFactor;
    
    %Refit:
    gcoefsX(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    gcoefsY(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))
    
    centers(:,i) = [ceil(gcoefsX(2,i)), ceil(gcoefsY(2,i))]; %center = [x y]
    sigmaX(:,i) = gcoefsX(3,i);
    sigmaY(:,i) = gcoefsY(3,i);
end

%Display every X fit:
if(0)
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsX(:,i),1:300));hold on; plot(mean(imageArrayC(:,:,i),1),'r'); hold off;      
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,4) == 0)       
        figure(i);
        plot(fg(gcoefsY(:,i),1:150));hold on; plot(mean(imageArrayC(:,:,i),2),'r'); hold off;      
    end
end
end

%%%%%Second moment:
xvector = [];
xvector = 1:length(imageArrayTC(1,:,1));
COMx = [];
for i=1:length(imageArrayTC(1,1,:))
        %Sum of Intensity*pixel location / sum of intensity
        COMx(i) = sum(mean(imageArrayTC(:,:,i),1).*xvector) / sum(mean(imageArrayTC(:,:,i),1));
end

%Second moment x direction
SMomX = []; sigmaXsm = [];
for i=1:length(imageArrayTC(1,1,:))
    SMomX(i) = sum(mean(imageArrayTC(:,:,i),1).*(xvector - COMx(i)).^2) / sum(mean(imageArrayTC(:,:,i),1));
end

sigmaXsm = sqrt(SMomX);

yvector = [];
yvector = 1:length(imageArrayTC(:,1,1));
COMy = [];
for i=1:length(imageArrayTC(1,1,:))
        %Sum of Intensity*pixel location / sum of intensity
        COMy(i) = sum(mean(imageArrayTC(:,:,i),2)'.*yvector) / sum(mean(imageArrayTC(:,:,i),2));
end

%Second moment y direction
SMomY = []; sigmaYsm = [];
for i=1:length(imageArrayTC(1,1,:))
    SMomY(i) = sum(mean(imageArrayTC(:,:,i),2)'.*(yvector - COMy(i)).^2) / sum(mean(imageArrayTC(:,:,i),2));
end

sigmaYsm = sqrt(SMomY);


%%%%%Atom numbers:
%Tight ROI array:
pixelCounts = [];
for i=1:length(imageArrayTC(1,1,:))
    pixelCounts(i) = sum(sum(imageArrayTC(:,:,i)));
end

%Pixel to Atom Correction:
pixelCounts(:) = varData(:,5)*0.42;


%Sort varData:
[sortedVarData,indexs] = sort(varData);
%indexs(:,1) is a vector of the sort.

pixelCountsSort = []; sigmaXSort = []; sigmaYSort = [];
TonTFsSort = []; imageArrayCSort = [];
sigmaRsmSort = []; TonTFsSort = [];

for i=1:length(sigmaX)
    sigmaXSort(i) = sigmaX(indexs(i));
    sigmaYSort(i) = sigmaY(indexs(i));
    pixelCountsSort(i) = pixelCounts(indexs(i));
    sigmaXsMomSort(i) = sigmaXsm(indexs(i));
    sigmaYsMomSort(i) = sigmaYsm(indexs(i));
    
    imageArrayCSort(:,:,i) = imageArrayC(:,:,indexs(i));
end

%Average over same BField data points:
j=1; runTotal = 0; magFields = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
imageArrayAvgs = [];
sMomentXStdDev = []; sMomentYStdDev = []; prev = sortedVarData(1,1);
for i=1:length(sortedVarData)
    curr = sortedVarData(i,1);
    
    if( curr == prev )
        runTotal = runTotal+1;
    else
        %hit next value
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        magFields(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        sMomentX(j) = mean(sigmaXsMomSort(i-runTotal:i));
        sMomentY(j) = mean(sigmaYsMomSort(i-runTotal:i));
        sMomentXStdDev(j) = std(sigmaXsMomSort(i-runTotal:i));
        sMomentYStdDev(j) = std(sigmaYsMomSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        runTotal = 0;
        j = j+1;
    end
    if( i == length(sortedVarData))
        %last run
        %disp('hit last run')
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        magFields(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        sMomentX(j) = mean(sigmaXsMomSort(i-runTotal:i));
        sMomentY(j) = mean(sigmaYsMomSort(i-runTotal:i));
        sMomentXStdDev(j) = std(sigmaXsMomSort(i-runTotal:i));
        sMomentYStdDev(j) = std(sigmaYsMomSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        
        runTotal = 0;
        j = j+1;
    end
    
    prev = curr;
end

%Fit functions to the averaged images:
gcoefsXa = []; gcoefsYa = []; gcoefsYaError = []; gcoefsXaError = [];
for i=1:length(imageArrayAvgs(1,1,:))
    [gcoefsXa(:,i),gcoefsXaError(:,:,i)] = gausFit1D(mean(imageArrayAvgs(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    [gcoefsYa(:,i),gcoefsYaError(:,:,i)] = gausFit1D(mean(imageArrayAvgs(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))       
end

widthsYavg = gcoefsYa(3,:);
widthsXavg = gcoefsXa(3,:);

widthsYavgError = []; widthsXavgError = [];
for i=1:length(widthsYavg)
    widthsYavgError(i) = widthsYavg(i) - gcoefsYaError(3,1,i);
    widthsXavgError(i) = widthsXavg(i) - gcoefsXaError(3,1,i);
end

%convert to real units:
%widthsX = widthsX.*pixelLength.*2; %*2 to make it not the radius
%widthsY = widthsY.*pixelLength.*2; 
%stdDevWidthsY = stdDevWidthsY.*pixelLength.*2; %full error on width
%stdDevWidthsX = stdDevWidthsX.*pixelLength.*2;

sMomentY = sMomentY.*pixelLength.*2; %*2 to make it not the radius
sMomentX = sMomentX.*pixelLength.*2; 
sMomentXStdDev = sMomentXStdDev.*pixelLength.*2; %full error on width
sMomentYStdDev = sMomentYStdDev.*pixelLength.*2;

figure(1);
h = errorbar(magFields,widthsY,stdDevWidthsY/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figname = [date '_' camera '_TransverseWidth_' num2str(atomnum) '_PixelNumber'];
figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\';
saveas(h,[figdirectory figname '.fig'],'fig');
saveas(h,[figdirectory figname '.png'],'png');
%line([832.2 832.2],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7]); 

figure(2);
errorbar(magFields,widthsX,stdDevWidthsX/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
figure(3);
errorbar(magFields,pixelNumbers,pixelNumbersStdDev/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
grid on;
if(0)
figure(4);
errorbar(magFields,sMomentY,sMomentYStdDev/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
figure(5);
errorbar(magFields,sMomentX,sMomentXStdDev/2,'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
figure(6);
plot(magFields,widthsX./(pixelNumbers.^0.25),'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
'Marker','o',...
'LineStyle','none',...
'Color',[0 0 1]);
end
        
figure(20)
errorbar(magFields,widthsYavg,widthsYavgError,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;
figure(21)
errorbar(magFields,widthsXavg,widthsXavgError,'MarkerSize',3,...
    'MarkerFaceColor',[0.600000023841858 0.600000023841858 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
grid on;




