%-------------------------------------------------------------------------%
%Load Images: This section is filename format dependent.

%Choose top or side camera:
topcam = 1;
sidecam = 0;

directory = 'C:\Data\140311_Imaging_Parameter_Check_2D_Trap_1usROI1_5usROI2_10usROI4_10usNewBeatROI5_1us20IsatROI3_Image_Pulse\';
datestring = '140311';
varstring = 'Isat';
varstring2 = 'HoldTime'; %if more than one relevant varstring or it changes throughout the log.
freqguess = 25; %2pi*x Hz guess

warning('OFF'); %#ok<WNOFF>
close all; %close all figures
%TotalImages = length(times)*shotsPerTime; %total number of images taken%


%Read in the log file:

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
                ftsImg = imread(fullname);
                imgArrayFresh(:,:,j) = ftsImg;
            else
                fail = 'fail on ROImax'
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

for i=1:length(fileLocList(:))
    
    if i<=101
        Isat = 135;
    elseif i<=151
        Isat = 135*5;
    elseif i<=161
        Isat = 135;
    else
        Isat = 135*10;
    end
    
   %disp(num2str(Isat));
   [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat); 
   imgArrayFresh(:,:,i) = atom2Image(:,:,1);

end



%-------------------------------------------------------------------------%
%Region of Interests:

%Region of Interests, use same for each image.
%Pulled from logfile:
ROIy = ROIbottom:ROItop;
ROIx = ROIleft:ROIright;

if(sidecam)
    ROIy = 276:323; %Top:Bot
    ROIx = 400:515;
end

if(topcam)
    ROIy = 38:88; %Top:Bot
    ROIx = 55:110;
end
ROIy = 10:170; %Top:Bot
ROIx = 10:200;

center = [(ROIx(length(ROIx))-ROIx(1))/2 + ROIx(1),(ROIy(length(ROIy))-ROIy(1))/2 + ROIy(1)]; %[x y]

%-------------------------------------------------------------------------%
%N_ROI Calc and sum over same Images:
 
%Sum over Isats (assumes contiguous within dataset):
sumImgArray = []; 
isatNums = []; %isatNums is [Isatval,number of]
currV = 0; currS = 0; %current variable (Isat) number of current set no
j=0;
for i=1:TotalImages
    if varData(i,1) == currV && varData(i,2) == currS
        isatNums(j,2) = isatNums(j,2) + 1;
        %Add another image to current element:
        %disp('Sum another image');
        for m=1:length(ROIy)
            for n=1:length(ROIx)
                sumImgArray(m,n,j) = sumImgArray(m,n,j) + imgArrayFresh(ROIy(m),ROIx(n),i);
            end
        end
                
    else
        j=j+1;        
        %disp('Add first image');
        currV = varData(i,1);
        currS = varData(i,2);
        %Add first element
        for m=1:length(ROIy)
            for n=1:length(ROIx)
                sumImgArray(m,n,j) = imgArrayFresh(ROIy(m),ROIx(n),i);
            end
        end
        %sumImgArray(:,:,j) = imgArrayFresh(ROIy,ROIx,i);
        isatNums(j,1) = currV;
        isatNums(j,2) = 1;
        isatNums(j,3) = currS;
        if currS == 1
            isatNums(j,4) = 1; %microsecond pulse length
        elseif currS == 2
            isatNums(j,4) = 5;
        elseif currS == 3
            isatNums(j,4) = 1;
        elseif currS == 4
            isatNums(j,4) = 10;
        else
            isatNums(j,4) = 10;
        end
    end
end

imgArray = []; %Image array reduced ie average over images


%Divide by the amount of images at that Isat to get average:

for i=1:length(sumImgArray(1,1,:))
    imgArray(:,:,i) = sumImgArray(:,:,i)/isatNums(i,2);
    %imgArray(:,:,i) = sumImgArray(:,:,i); %To not average, just use sum.
end

%Noise correction from small region beside cloud:

NoiseXL = 5:10; NoiseXR = 155:160; NoiseXT = 5:10;
for i=1:length(imgArray(1,1,:))
    avgNoiseL = mean(mean(imgArray(:,NoiseXL,i)));
    avgNoiseR = mean(mean(imgArray(:,NoiseXR,i)));
    avgNoiseT = mean(mean(imgArray(NoiseXT,:,i)));
    %avgNoiseL
    %avgNoiseR
    avgNoise = (avgNoiseL+avgNoiseR)/2;
    imgArray(:,:,i) = imgArray(:,:,i) - avgNoiseT;
end

%Order each dataset:

[isatSorted,indexs] = sortrows(isatNums,[4 3 1]);
sortArray = [];

for i=1:length(imgArray(1,1,:))
    sortArray(:,:,i) = imgArray(:,:,indexs(i));
end

imgArray = sortArray;

%Region of Interest Number:

NROISum = [];
for i=1:length(imgArray(1,1,:))
    NROISum(i) = sum(sum(imgArray(:,:,i)));
end

%Profiles:

profArray = [];
for i=1:length(imgArray(1,1,:))
    %profArray(:,i) = sortArray(:,length(ROIx)/2,i);%Single slice in mid
    %profArray(:,i) = sum(imgArray(1:length(ROIy),:,i));%Horz. slice
    profArray(:,i) = sum(imgArray(:,1:length(ROIx),i));%Vert. slice
end

%Guassian fit to profiles:
gfitcoefs = [];

for i=[7,28]
%Fit:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4)); %function to fit with
p0 = [700 85 20 5];
lb = [400 40 10 1];
ub = [1000 100 50 30];

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
xs = 1:length(profArray(:,i));
coefs = lsqcurvefit(fg,p0,xs(:),profArray(:,i),lb,ub,curvefitoptions);
%[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs(:),profArray(:,7),lb,ub,curvefitoptions);

gfitcoefs(:,i) = coefs;

end


%Plot Fit:
%plot(xs,fg(coefs,xs),'r');


%Side by Side plots:
%for i=1:length(imgArray(1,1,:))
if(0) %Skip this
for i=[7,28]
    
    titular = ['N = ' num2str(NROISum(i),3) ...
        ', Pulse = ' num2str(isatSorted(i,4)) 'us, Isat*' ...
            num2str(isatSorted(i,1)) ', Average of ' num2str(isatSorted(i,2)) ' shots'];

    figure('Position',[400, 150, 600, 800]);
    subplot(2,1,1);image(imgArray(:,:,i),'CDataMapping','scaled'); caxis ([-1 15]);
    title(titular); 
    subplot(2,1,2); hold on; plot(profArray(:,i)); plot(xs,fg(gfitcoefs(:,i)',xs),'r'); hold off;  
       

end
end

%Special Gaussian conversion plot:


figure('Position',[200, 150, 1600, 800]); %[left,bot,width,height]
%7
i=7;   
titular7 = ['N = ' num2str(NROISum(i),3) ...
    ', Pulse = ' num2str(isatSorted(i,4)) 'us, Isat*' ...
    num2str(isatSorted(i,1)) ', Average of ' num2str(isatSorted(i,2)) ' shots'];
subplot(2,3,1);image(imgArray(:,:,i),'CDataMapping','scaled'); caxis ([-1 15]);
title(titular7);
subplot(2,3,4); hold on; plot(profArray(:,i)); plot(xs,fg(gfitcoefs(:,i)',xs),'r'); hold off;
title('Red: Gaussian fit to profile of above cloud.');
%28
i=28;
titular28 = ['N = ' num2str(NROISum(i),3) ...
    ', Pulse = ' num2str(isatSorted(i,4)) 'us, Isat*' ...
    num2str(isatSorted(i,1)) ', Average of ' num2str(isatSorted(i,2)) ' shots'];
subplot(2,3,3);image(imgArray(:,:,i),'CDataMapping','scaled'); caxis ([-1 15]);
title(titular28);
subplot(2,3,6); hold on; plot(profArray(:,i)); plot(xs,fg(gfitcoefs(:,i)',xs),'g'); hold off;
title('Green: Gaussian fit to profile of above cloud.');
subplot(2,3,5); %bottom middle plot
hold on; 
plot(xs,fg(gfitcoefs(:,7)',xs),'r'); 
plot(xs,fg(gfitcoefs(:,28)',xs),'g'); 
%Scale gaussian:
scalecoefs(1) = gfitcoefs(1,28)*1.5; scalecoefs(2) = gfitcoefs(2,7); scalecoefs(3) = gfitcoefs(3,28); scalecoefs(4) = gfitcoefs(4,28);
plot(xs,fg(scalecoefs,xs),'black');
title('Black: Attempted scaling of green to red. Multiply coef by 1.5.');
hold off;


%Spectrum Binning of column density values -------------------------------%
if(0)
imageToCalc = 7;
binMaxSize = length(ROIx)*length(ROIy)*1.2;
spectrumBin = [];

for i=1:binMaxSize
    spectrumBin(i) = 0;
end


maxValue = max(max(imgArray(:,:,imageToCalc))); %maxValue of image
spectrumBin(end) = maxValue;

for i=1:length(ROIx)
    for j=1:length(ROIy)
        toPut = imgArray(j,i,imageToCalc);
        %Removes all less than zero elements...
        if toPut < 0 
            toPut = 0;
        end
        relativeLoc = ceil(toPut/maxValue * binMaxSize);
        
        %Only put in if that location is empty (non-filled)
        if(relativeLoc > 0)
            if(spectrumBin(relativeLoc) ~= 0)
                spectrumBin(relativeLoc) = toPut;
            end
        end
    end
end

%Empty out the zeros in the array:
spectrumBin(spectrumBin == 0) = [];
end


%Column density spectrum conversion---------------------------------------%

%From 7 to 28:
convFrom = 7;
convTo = 28;












%-------------------------------------------------------------------------%
%Fancy images plot%

collated = [];
textLocBR = [];

while j*length(ROIx) <= 500
    j = j+1;
end

intWidths = j;
k=0;
total = 1;
i=0; inc = 1;
while total <= length(imgArray(1,1,:))
    
        collated(k*length(ROIy)+mod(1,k):(k*length(ROIy)+length(ROIy)-1+mod(1,k)),...
        i*length(ROIx)+mod(1,i):(i*length(ROIx)+length(ROIx)-1+mod(1,i)),1)...
        = imgArray(1:length(ROIy),1:length(ROIx),total);

        textLocBR(total,1) = (k*length(ROIy)+length(ROIy)-1+mod(1,k)) - 10; %Textpos y(40x40 box)
        textLocBR(total,2) = (i*length(ROIx)+length(ROIx)-1+mod(1,i)) - 40;  %Textpos x(40x40 box)
    
    
        i = i+1;
        total = total+1;

    
        %if mod(i,ceil(length(imgArray(1,1,:))/intWidths)+1) == 0
        %    disp('downshift');
        %    k = k+1;         
        %   i = 1;
        %end
         
        if total <= length(imgArray(1,1,:))
            if isatSorted(total-1,3) ~= isatSorted(total,3)
                %disp('downshift');
                k = k+1;
                i = 0;
            end
        end
  
end

figure(98);
%imagesc(collated);
image(collated,'CDataMapping','scaled'); caxis ([-1 10]);

%Apply text:
for i=1:length(textLocBR(:,1))
   
    %t = text(textLocBR(i,2),textLocBR(i,1),[num2str(isatSorted(i,4),2) 'us'], ...
         %'FontSize',6);  
end

hold off;


%-------------------------------------------------------------------------%
%Summed Image:%
sumImages = [];

    %Whole Image Sum (Slow)
    %for j=1:length(imgArray(:,1,1))
    %    for k=1:length(imgArray(1,:,1))
    %sumImages(j,k) = sum(imgArray(j,k,:));
    %    end
    %end
    
    for j=1:length(ROIy)
        for k=1:length(ROIx)
    sumImages(j,k) = sum(imgArray(j,k,:));
        end
    end
    
figure(100);
scatter(1:length(imgArray(1,1,:)),NROISum)
    
figure(99);
image(sumImages,'CDataMapping','scaled');
hold off;

%image(imgArray(:,:,20),'CDataMapping','scaled'); caxis ([-1 20]); isatSorted(20,:)
