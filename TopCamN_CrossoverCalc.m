function [logSigmaY,logSigmaX,logNROISum] = TopCamN_CrossoverCalc(directoryI,datestringI,binsI)
%-------------------------------------------------------------------------%
%Initialization Area:

%Choose top or side camera:
topcam = 1;
sidecam = 0;

if(topcam)
    fittingExp = 0; %only for non-binned
    fittingGaus = 0; %Fitting gaussian method
    secondMoment = 1; %Second moment method
    dispNumber = 6; %display every xth image
    bin = 1; %binning on or off
    %bins = 16; %number of bins
    bins = binsI;
    centering = 1; %centering bins on or off
    leveling = 1; %noise level correction
    figures = 0; %Suppress figures
    resetFigures = 0; %wipe all existing figures
elseif(sidecam)
    fittingExp = 0; %only for non-binned
    fittingGaus = 0; %Fitting gaussian method
    secondMoment = 0; %Second moment method
    dispNumber = 6; %display every xth image
    bin = 0; %binning on or off
    %bins = 16; %number of bins
    bins = binsI;
    centering = 0; %centering bins on or off
    leveling = 0; %noise level correction
    figures = 1; %Suppress figures
    resetFigures = 1; %wipe all existing figures
end

%directory = 'C:\Data\140318_880G_Crossover_Measurement_10us_Imagepulse_1_Isat\';
directory = directoryI;
%datestring = '140318';
datestring = datestringI;
varstring = 'EvapVoltage';
varstring2 = 'HoldTime'; %if more than one relevant varstring or it changes throughout the log.
Isat = 135*10; %Isat value, 135 for 1us
pi = 3.14159;

warning('OFF'); %#ok<WNOFF>
if(resetFigures)
close all; %close all figures
end
%TotalImages = length(times)*shotsPerTime; %total number of images taken%



%-------------------------------------------------------------------------%
%Load Images: This section is filename format dependent.

%Read in the log file:
disp('Reading in log file...');

logfilename = [directory datestring '_log_camera.txt'];
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

imgArrayFresh = []; 
fileLocList = cell(TotalImages,1);
prev = 0;
j=0; k=0;
dataset = 1;

%for loop to get information out of the logfile and load it into arrays:
for i=1:(length(C)-1)
    
    if i>1
        prev = C{i-1};
    end
    curr = C{i};
    next = C{i+1};
    %ROImax = C{i+20};
    ROImax = '1';
    
    %load sidecam images
    if(sidecam)
        if strcmp(curr,'MeasNr')
            if str2num(ROImax) ~= 0
                j = j+1;
                filename = prev;
                picno = ['_' filename(end-2:end)];
                fullname = [directory filename picno '.fts'];
                %fullname = [directory filename '_001' '.fts'];
                %ftsImg = imread(fullname);
                fileLocList{j} = fullname;
                %imgArrayFresh(:,:,j) = ftsImg;
            else
                disp('fail on ROImax');
            end
        end
    end
    
    %load topcam images
    if(topcam)
        if strcmp(curr,'#WinView#')
            j = j+1;
            filename = next;
            fileLocList{j} = [directory filename];
            %[beamImage,atom1Image,atom2Image] = PullSPE(directory,filename);
            %imgAtom2Fresh(:,:,j) = atom2Image(:,:,1); %state 2 images in atom2Image
        end
    end
    
    %data set identifier
    dataset = 1;
    if strcmp(curr,'N_ROI')
        dataset = 1;
    elseif strcmp(curr,'N_ROI2')
        dataset = 2;
    elseif strcmp(curr,'N_ROI3')
        dataset = 3;
    elseif strcmp(curr,'N_ROI4')
        dataset = 4;
    elseif strcmp(curr,'N_ROI5')
        dataset = 5;
    end
    
    if strcmp(curr,varstring) || strcmp(curr,varstring2)
        k=k+1;
        %logholdtime(ceil(i/(length(C)/TotalImages))) = str2num(next);
        varData(k,2) = dataset;
        varData(k,1) = str2num(next);
    end
    
    %Also pull ROI out:
    if strcmp(curr,'ROItop')
        ROItop = str2num(next);
    elseif strcmp(curr,'ROIbottom')
        ROIbottom = str2num(next);
    elseif strcmp(curr,'ROIleft')
        ROIleft = str2num(next);
    elseif strcmp(curr,'ROIright')
        ROIright = str2num(next);
    end
       

end  
%-------------------------------------------------------------------------%
%From list of .SPE filenames, call PullSPE.m which returns the images
%converted to Atom Number already with all corrections. 

if(topcam)
    disp('Processing SPE Images...');
    
    for i=1:length(fileLocList(:))
        
        %disp(num2str(Isat));
        [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,0);
        imgArrayFresh(:,:,i) = atom2Image(:,:,1);
        
    end
end

%From list of .FTS filenames, call PullFTS.m which returns the images
%with all corrections. 


if(sidecam)
    ftsImage = [];
    disp('Processing FTS Images...');
    
    for i=1:length(fileLocList(:))
        
        [ftsImage] = PullFTS(fileLocList{i});
        imgArrayFresh(:,:,i) = ftsImage;
        
    end
    
end



%-------------------------------------------------------------------------%
%Region of Interests:

%Region of Interests, use same for each image.
%Pulled from logfile:
%ROIy = ROIbottom:ROItop;
%ROIx = ROIleft:ROIright;

if(sidecam)
    ROIy = 500:650; %Top:Bot
    ROIx = 420:720;
end

if(topcam)
%Wide region
ROIy = 10:170; %Top:Bot
ROIx = 10:200;

if(secondMoment)
%Tight region
ROIy = 30:140; %Top:Bot
ROIx = 30:170;
end
end


%ROIy = 1:(length(imgArrayFresh(:,1,1))-20);
%ROIx = 1:(length(imgArrayFresh(1,:,1))-20);


%center = [(ROIx(length(ROIx))-ROIx(1))/2 + ROIx(1),(ROIy(length(ROIy))-ROIy(1))/2 + ROIy(1)]; %[x y]

disp('Established ROIs...');

%Atom Number, Profiles, Binning-------------------------------------------%

imgArray = []; %Image array reduced ie average over images

if(bin ~= 1)
    %Define imgArray:
    for i=1:length(imgArrayFresh(1,1,:))
        imgArray(:,:,i) = imgArrayFresh(ROIy,ROIx,i);
    end
    
        
    %Atom Number:
    
    %Atom number for ROISum on each individual image.
    NROISum = [];
    for i=1:length(imgArray(1,1,:))
        NROISum(i) = sum(sum(imgArray(:,:,i)));
    end
    maxNumber = max(NROISum) + 1;
end


%----Pre-fit Gaussians to all data for centering etc-----%
%This section goes through each individual image (before binning)
%and fits a gaussian to find the peak location; and hence the center
%of each cloud. It then uses this information to shift the binning
%to the centre of the first cloud in that bin.

disp('Fitting Gaussians to All...');

imgYFs = []; %Array of the Y profiles
imgXFs = [];
gfitcoefsYF = []; gfitcoefsXF = []; %Arrays of fit co-efficients
xsX = []; xsY = [];
profArrayY = [];

%Y (vertical) Calculations:
%Sum images along y direction:
for i=1:length(imgArrayFresh(1,1,:))
    for j=1:(length(imgArrayFresh(:,1,1)))
        imgYFs(j,i) = sum(imgArrayFresh(j,:,i)); %Sum each horizontal 1D array
    end
end

profArrayY = imgYFs;

%Fit gaussian and find index of the peak
for i=1:length(imgArrayFresh(1,1,:))
    %Fit:
    fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4)); %function to fit with
    if(topcam)
        p0 = [700 85 20 5];
        lb = [50 50 10 -30];
        ub = [1400 110 80 30];
    elseif(sidecam)
        p0 = [160000 575 20 5];
        lb = [30000 520 0.01 -5000];
        ub = [300000 650 120 5000];
    end
    
    
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    xsY = 1:length(imgYFs(:,i));
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xsY(:),profArrayY(:,i),lb,ub,curvefitoptions);
    %[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArray(:,7),lb,ub,curvefitoptions);
    
    gfitcoefsYF(:,i) = coefs;
    %gfitciY(:,:,i) = nlparci(coefs,r,'jacobian',J);
    
end

%plot(fg(gfitcoefsYF(:,5),xsY),'r'); hold on; plot(imgYFs(:,5),'b');
profArrayX = [];

%X (horizontal) Calculations:
%Sum images along x direction:
for i=1:length(imgArrayFresh(1,1,:))
    for j=1:(length(imgArrayFresh(1,:,1)))
        imgXFs(j,i) = sum(imgArrayFresh(:,j,i)); %Sum each horizontal 1D array
    end
end

profArrayX = imgXFs; %The X profiles

for i=1:length(imgArrayFresh(1,1,:))
    %Fit:
    fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4)); %function to fit with
    if(topcam)
        p0 = [700 110 20 5];
        lb = [50 90 10 -200];
        ub = [1400 130 80 50];
    elseif(sidecam)
        p0 = [25000 575 20 5];
        lb = [1000 520 10 -100];
        ub = [30000 650 120 100];
    end
    
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    xsX = 1:length(profArrayX(:,i));
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xsX(:),profArrayX(:,i),lb,ub,curvefitoptions);
    %[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArray(:,7),lb,ub,curvefitoptions);
    
    gfitcoefsXF(:,i) = coefs;
    %gfitciX(:,:,i) = nlparci(coefs,r,'jacobian',J);
    
end

indexY = []; indexX = [];
for i=1:length(imgArrayFresh(1,1,:))
    [valueY, indexY(i)] = max(fg(gfitcoefsYF(:,i),xsY));
    [valueX, indexX(i)] = max(fg(gfitcoefsXF(:,i),xsX));
end

tROIy = []; tROIx = [];
%Really tight dynamic ROI's for 2nd moment calc
for i=1:length(imgArrayFresh(1,1,:))
    tROIy(1,i) = indexY(i) - round(gfitcoefsYF(3,i)*1); %1 times the std deviation
    tROIy(2,i) = indexY(i) + round(gfitcoefsYF(3,i)*1);
    tROIx(1,i) = indexX(i) - round(gfitcoefsXF(3,i)*1);
    tROIx(2,i) = indexX(i) + round(gfitcoefsXF(3,i)*1);
end


if(leveling)
    disp('Leveling...');
    %Re-level the images so that the background is around zero:
    %ie extra background noise subtraction...
    shift = [];
    %Noise area:
    NROIy = 160:190;
    NROIx = 40:160;
    
    for i=1:length(imgArrayFresh(1,1,:))
        %shift = -(min(gfitcoefsXF(4,i),gfitcoefsYF(4,i)));
        shift(i) = mean(mean(imgArrayFresh(NROIy,NROIx,i)));
        for j=1:length(imgArrayFresh(:,1,1))
            for k=1:length(imgArrayFresh(1,:,1))
                imgArrayFresh(j,k,i) = imgArrayFresh(j,k,i) - shift(i);
            end
        end
    end
end


%Atom number for ROISum on each individual image.
NROISumImage = [];
for i=1:length(imgArrayFresh(1,1,:))
    NROISumImage(i) = sum(sum(imgArrayFresh(ROIy,ROIx,i)));
end

maxNumber = max(NROISumImage) + 1;
    
if(bin)
    disp('Binning Data...');    
    
   
    bROIy = []; bROIx = []; %Tight regions for each bin for 2nd moment calc.
    
    trueCenterY = []; trueCenterX = [];
    binArray = [];
    binCounts = [];
    binNumber = []; %Which bin the image ends up in.
    %bins = 24;
    binned = 0;
    
    for i=1:bins
        for m=1:length(ROIy)
            for n=1:length(ROIx)
                binArray(m,n,bins) = 0;
            end
        end
    end
    
    for i=1:bins
        binCounts(i) = 0;
    end
    

    j=1;
    for i=1:length(imgArrayFresh(1,1,:))

        while(binned == 0)
            if NROISumImage(i) <= (j/bins)*maxNumber               
                %Calculate Center offset relative to first image in
                %binned set:
                if(centering)
                    if(binCounts(j) == 0)
                        %first image in this bin
                        trueCenterY(j) = indexY(i);
                        trueCenterX(j) = indexX(i);
                        diffY = 0;
                        diffX = 0;
                        bROIy(1,j) = tROIy(1,i);
                        bROIy(2,j) = tROIy(2,i);
                        bROIx(1,j) = tROIx(1,i);
                        bROIx(2,j) = tROIx(2,i);
                    else
                        diffY = trueCenterY(j) - indexY(i);
                        diffX = trueCenterX(j) - indexX(i);
                    end
                
                if(diffX > 10 || diffY > 10 || diffX < -10 || diffY < -10)
                    %diff is very large(!)
                    %ie if it's greater than 10 pixels movement ignore. 
                    disp('Warning: Diff is very large.');
                    %diffX = 0;
                    %diffY = 0;
                    disp([num2str(diffY) '=diffy ' num2str(diffX) '=diffx ']);
                end
                end
                
                %disp([num2str(diffY) '=diffy ' num2str(diffX) '=diffx ']);
                %Add images elementwise
                for m=1:length(ROIy)
                    for n=1:length(ROIx)
                        if(centering)
                            %making sure real and positive integers
                            %and within the larger pullSPE regions
                            if(ROIy(m) + diffY <= 0)
                                diffY = 0;
                            end
                            if(ROIx(n) + diffX <= 0)
                                diffX = 0;
                            end
                            if(ROIx(n) + diffX >= 206)
                                diffX = 0;
                            end
                            if(ROIy(m) + diffY >= 250)
                                diffY = 0;
                            end
                            
                        binArray(m,n,j) = binArray(m,n,j)...
                            + imgArrayFresh(ROIy(m) + diffY,ROIx(n) + diffX,i);
                        else
                        binArray(m,n,j) = binArray(m,n,j)...
                            + imgArrayFresh(ROIy(m),ROIx(n),i);    
                        end
                        
                    end
                end
                %binArray(:,:,j) = binArray(:,:,j) + imgArrayFresh(:,:,i);
                binCounts(j) = binCounts(j) + 1;
                binNumber(i) = j;
                
                binned = 1;
            end
            j = j+1;
        end
        
        binned = 0;
        j=1;  
        
    end
    
    
    for i=1:bins
        if binCounts(i) == 0
            disp('Bin with Zero elements found!');
            imgArray(:,:,i) = binArray(:,:,i);
        else
            imgArray(:,:,i) = binArray(:,:,i)./binCounts(i);
        %imgArray(:,:,i) = binArray(:,:,i);
        end
    end
    
    
    
    %Atom number for ROISum on each bin image.
    NROISum = [];
    for i=1:length(imgArray(1,1,:))
        NROISum(i) = sum(sum(imgArray(:,:,i)));
    end
end


%Second moment calc: -----------------------------------------------------%

%Constrain ROI's:
%if(secondMoment)
    %ROIy = min(tROIy(1,:)):max(tROIy(2,:));
    %ROIx = min(tROIx(1,:)):max(tROIx(2,:));
%end



%Centre-of-mass calc:
imgYs = []; %Array of the Y profiles (heights) to find COM movement on.
COMy = []; %Array of Image centre of masses of the y profiles.
imgXs = [];
COMx = [];
imgYm = [];
imgXm = [];

%VERTICAL CALCS (Y)%%%%%%%%%%%%%%%%%%
%Sum images along x direction:
for i=1:length(imgArray(1,1,:))
    for j=1:length(ROIy)
            imgYs(j,i) = sum(imgArray(j,:,i)); %Sum each horizontal 1D array
            imgYm(j,i) = mean(imgArray(j,:,i)); %Avg each horizontal 1D array
    end
end

if(secondMoment)
disp('Calculating Moments...');
%Find COM in y direction:
yvector = 1:length(ROIy);

for i=1:length(imgArray(1,1,:))
        %Sum of Intensity*pixel location / sum of intensity
        COMy(i) = sum(imgYs(:,i).*yvector(:)) / sum(imgYs(:,i));
end


    %Second moment y direction
    SMomY = [];
    for i=1:length(imgArray(1,1,:))
        SMomY(i) = sum(imgYs(:,i).*(yvector(:) - COMy(i)).^2) / sum(imgYs(:,i));
    end
    
    sigmaY = sqrt(SMomY);
end

%HORIZONTAL CALCS (X)%%%%%%%%%%%%%%%%%%
%Sum images along x direction:
for i=1:length(imgArray(1,1,:))
    for j=1:length(ROIx)
            imgXs(j,i) = sum(imgArray(:,j,i)); %Sum each horizontal 1D array
            imgXm(j,i) = mean(imgArray(:,j,i)); %Avg each horizontal 1D array
    end
end

if(secondMoment)
%Find COM in x direction:
xvector = 1:length(ROIx);

for i=1:length(imgArray(1,1,:))
        %Sum of Intensity*pixel location / sum of intensity
        COMx(i) = sum(imgXs(:,i).*xvector(:)) / sum(imgXs(:,i));
end


    %Second moment x direction
    SMomX = [];
    for i=1:length(imgArray(1,1,:))
        SMomX(i) = sum(imgXs(:,i).*(xvector(:) - COMx(i)).^2) / sum(imgXs(:,i));
    end
    
    sigmaX = sqrt(SMomX);
end


%Fitting: ----------------------------------------------------------------%

disp('Fitting functions...');

if(fittingGaus)
    gfitcoefsY = []; gfitcoefsX = [];
    gfitciY = []; gfitciX = [];
    
    %Fit to Y or X:
    profArrayYm = imgYm; %First Y profiles
    sigmaY = []; sigmaX = [];
    
    for i=1:length(imgArray(1,1,:))
        %Fit:
        fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4)); %function to fit with
        p0 = [4 85 20 5];
        lb = [0.1 40 10 -30];
        ub = [30 180 50 30];
        
        %fg = @(p,x)(1./(sqrt(2*pi)*p(1)).*exp((-1).*((x-p(2)).^2) ./ (2.*p(1).^2)) + p(3)); %function to fit with
        %p0 = [20 700 5];
        %lb = [10 50 -30];
        %ub = [50 1000 30];
        
        curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
        xs = 1:length(profArrayYm(:,i));
        [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArrayYm(:,i),lb,ub,curvefitoptions);
        %[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArray(:,7),lb,ub,curvefitoptions);
        
        gfitcoefsY(:,i) = coefs;
        gfitciY(:,:,i) = nlparci(coefs,r,'jacobian',J);
        sigmaY = gfitcoefsY(3,:); %Standard deviation
                   
    end
    
    profArrayXm = imgXm; %Then X profiles
    
    for i=1:length(imgArray(1,1,:))
        %Fit:
        fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4)); %function to fit with
        p0 = [4 85 20 5];
        lb = [0.1 40 10 -100];
        ub = [30 180 50 50];
        
        %fg = @(p,x)(1./(sqrt(2*pi)*p(1)).*exp((-1).*((x-p(2)).^2) ./ (2.*p(1).^2)) + p(3)); %function to fit with
        %p0 = [20 700 5];
        %lb = [10 50 -30];
        %ub = [50 1000 30];
               
        curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
        xs = 1:length(profArrayXm(:,i));
        [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArrayXm(:,i),lb,ub,curvefitoptions);
        %[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArray(:,7),lb,ub,curvefitoptions);
        
        gfitcoefsX(:,i) = coefs;
        gfitciX(:,:,i) = nlparci(coefs,r,'jacobian',J);
        sigmaX = gfitcoefsX(3,:); %Standard deviation
                   
    end
    
    
end

if(fittingExp)
    fg = @(p,x)(p(1).*x.^(p(2))); %function to fit with
    p0 = [6 1/4];
    lb = [1 1/10];
    ub = [8 1/2];
    
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    xs = 1:length(profArray(:,i));
    
    sepIndex = 1;
    %Order NROISum:
    [SortedNROISum,indexs] = sortrows(NROISum');
    %NROISum = SortedNROISum';
    
    %Order sigmax, sigmay:
    for i=1:length(sigmaX)
        sigmaXSort(i) = sigmaX(indexs(i));
        sigmaYSort(i) = sigmaY(indexs(i));
    end
    
    %Split the data up to ostensibly 2D and 3D
    for i=1:length(SortedNROISum)
        if SortedNROISum(i) > 25000
            display('SplitData Success');
            sepIndex = i;
            break;
        end
    end
    
    firstHalfData = 20:sepIndex;
    secondHalfData = sepIndex+1:length(NROISum);
    
    coefsX1 = lsqcurvefit(fg,p0,SortedNROISum(firstHalfData)',sigmaXSort(firstHalfData),lb,ub,curvefitoptions);
    coefsY1 = lsqcurvefit(fg,p0,SortedNROISum(firstHalfData)',sigmaYSort(firstHalfData),lb,ub,curvefitoptions);
    coefsX2 = lsqcurvefit(fg,p0,SortedNROISum(secondHalfData)',sigmaXSort(secondHalfData),lb,ub,curvefitoptions);
    coefsY2 = lsqcurvefit(fg,p0,SortedNROISum(secondHalfData)',sigmaYSort(secondHalfData),lb,ub,curvefitoptions);
    
    fullxs = 1:10:max(NROISum); %Full range of ROI Sum to plot against
end

if(~bin && topcam)
    SortedNROISum = []; sigmaXSort = []; sigmaYSort = [];
   %Sort the arrays:
    [SortedNROISum,indexs] = sortrows(NROISum');
    %NROISum = SortedNROISum';
    
    %Order sigmax, sigmay:
    for i=1:length(sigmaX)
        sigmaXSort(i) = sigmaX(indexs(i));
        sigmaYSort(i) = sigmaY(indexs(i));
    end
    
    sigmaY = sigmaYSort;
    sigmaX = sigmaXSort;
    NROISum = SortedNROISum';
    
end

if(sidecam)
    %For the sidecam, we use the initial image fits for the sigma data:
    sigmaX = gfitcoefsXF(3,:); %Standard deviation
    sigmaY = gfitcoefsYF(3,:); 
    NROISum = NROISumImage;
end



%logs of data:
logNROISum = log10(NROISum);
logSigmaY = log10(sigmaY);
%logSigmaY = sigmaY;
logSigmaX = log10(sigmaX);

%Plotting: ---------------------------------------------------------------%

disp('Plotting...');
%if figures = 0 suppress output
if(figures)

%Side by Side plots:
%for i=1:length(imgArray(1,1,:))

if(sidecam)
    
    
    for i=1:length(imgArray(1,1,:))
        if(mod(i,dispNumber) == 0) %Only disp every one...
            
            figure('Position',[400, 150, 1500, 800]);
            subplot(2,2,1);imagesc(imgArrayFresh(:,:,i));%,'CDataMapping','scaled'); caxis ([-1 15]);
            subplot(2,2,2); hold on; plot(imgYFs(:,i));
            plot(xsY,fg(gfitcoefsYF(:,i)',xsY),'r');
            subplot(2,2,3); hold on; plot(imgXFs(:,i));
            plot(xsX,fg(gfitcoefsXF(:,i)',xsX),'r');
            hold off;
            
        end
    end
    
    
    
      
    
else

if(fittingGaus)
    %figure('Position',[400, 150, 600, 800]);
    for i=1:length(imgArray(1,1,:))
        if(mod(i,dispNumber) == 0) %Only disp every one...
        
        figure('Position',[400, 150, 600, 800]);
        subplot(2,1,1);image(imgArray(:,:,i),'CDataMapping','scaled'); caxis ([-1 15]);
        %title(titular);
        subplot(2,1,2); hold on; plot(profArray(:,i)); 
        plot(xs,fg(gfitcoefsX(:,i)',xs),'r'); 
        hold off;
        
        end 
    end

end

%Summed Image:%
sumImages = [];

for j=1:length(ROIy)
    for k=1:length(ROIx)
        sumImages(j,k) = sum(imgArray(j,k,:));
    end
end


figure(101);
image(sumImages,'CDataMapping','scaled');
title('Sum over all images in ROI');

figure(100);
scatter(1:length(imgArray(1,1,:)),NROISum);
title('Image Number vs Atom Number');


%figure(97);
%scatter(log(NROISum),log(SMom));

%figure(96);
%scatter(NROISum,log(SMom));



if(bin ~= 1)    
figure(99);
scatter(varData(:,1),NROISum);
title('Evaporation mV vs Atom Number');

figure(95);
hold on;
scatter(NROISum,SMomY,'b');
scatter(NROISum,SMomX,'r');
hold off;
title('Atom Number vs Second Moment. Red=X Blue=Y');
  
figure(94);
hold on;
scatter(varData(:,1),sigmaY,'b');
scatter(varData(:,1),sigmaX,'r');
hold off;
title('Evaporation mV vs Sigma. Red=X Blue=Y');
end

if(fittingExp)
figure(93);
hold on;
scatter(NROISum,sigmaY,15,[0 0 1]); %blue [0 0 1]
%scatter(NROISum,sigmaX,15,[0.5 0 0]);
%plot(fullxs,fg(coefsX1,fullxs),'Color',[0.7 0 0]);
%plot(fullxs,fg(coefsX2,fullxs),'Color',[1 0 0]); %red [1 0 0]
plot(fullxs,fg(coefsY1,fullxs),'Color',[0 0 0.3]);
plot(fullxs,fg(coefsY2,fullxs),'Color',[0 0 0.6]);
hold off;
title('Atom Number vs Sigma. Red=X Blue=Y');
elseif(fittingGaus)
    errorSigmaXlb = []; errorSigmaXub = [];
    errorSigmaYlb = []; errorSigmaYub = [];
    
    %Calculate error bars for gaussian fit:
    for i=1:length(sigmaX)
        errorSigmaXlb(i) = sigmaX(i)-gfitciX(3,1,i);
        errorSigmaXub(i) = gfitciX(3,2,i) - sigmaX(i);
    end
    
    for i=1:length(sigmaY)
        errorSigmaYlb(i) = sigmaY(i)-gfitciY(3,1,i);
        errorSigmaYub(i) = gfitciY(3,2,i) - sigmaY(i);
    end
    
figure(111);
errorbarxy(NROISum,sigmaY,zeros(20),errorSigmaYlb,zeros(20),errorSigmaYub,'.');
title('Atom Number vs Sigma Y'); 
figure(112);
errorbarxy(NROISum,sigmaX,zeros(20),errorSigmaXlb,zeros(20),errorSigmaXub,'.');
title('Atom Number vs Sigma X'); 

else
 figure(93);
hold on;
scatter(NROISum,sigmaY,'b','o');
scatter(NROISum,sigmaX,'r','o');
hold off;
title('Atom Number vs Sigma. Red=X Blue=Y');

end


end
%end if(figures)
end
%r =
%rscan(imgArray(:,:,3),'xavg',length(imgArray(1,:,3))/2,'yavg',length(imgArray(:,1,3))/2);


%end function
end