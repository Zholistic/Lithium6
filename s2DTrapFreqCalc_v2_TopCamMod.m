%2D Trap Side Imaging Frequency Calculation Script
%Tyson Peppler jan-2014
%
% March 2014 mod to handle TopCam .SPE files.
%
%-------------------------------------------------------------------------%
%Load Images: This section is filename format dependent.

%Choose top or side camera:
topcam = 0;
sidecam = 1; raw = 1;

directory = 'C:\Data\140414_SIDECAM_MeasureTrapFrequency_2ndTrap_P1W_AOMDriver9V4\';
datestring = '140414';
varstring = 'HoldTime';
varstring2 = 'Holdtime'; %if more than one relevant varstring or it changes throughout the log.
freqguess = 25; %2pi*x Hz guess
Isat = 135*10;

warning('OFF'); %#ok<WNOFF>
close all;
%TotalImages = length(times)*shotsPerTime; %total number of images taken%


%Read in the log file:

logfilename = [directory datestring '_log_camera.txt'];
fid = fopen(logfilename,'rt');
C = textscan(fid, '%s', 'Delimiter','\t'); %tokenize into tab seperated tokens
C = C{1};
fclose(fid);

%Iterate through strings in the log file array 1-d C{:} to find the
%information we want:

logholdtime = [];

TotalImages = 0;

%for loop to find the total number of images
for i=1:(length(C)-1)

    curr = C{i};
    next = C{i+1};
    
    if strcmp(curr,'MeasNr')
        TotalImages = TotalImages + 1;
    end
              
    
end

fileLocList = cell(TotalImages,1);
imgArray = []; 
fileList = [];
prev = 0;
j=0; k=0;

%for loop to get information out of the logfile and load it into arrays:
for i=1:(length(C)-1)
    
    if i>1
        prev = C{i-1};
    end
    curr = C{i};
    next = C{i+1};
    %ROImax = C{i+20};
    ROImax = '1';
    
    if(0)
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
                imgArray(:,:,j) = ftsImg;
            else
                fail = 'fail on ROImax'
            end
        end
    end
    end
    
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
            %fileList(j) = cellstr(filename);
            fileloc = [directory filename];
            [beamImage,atom1Image,atom2Image] = PullSPE(fileloc,Isat,0); %OD = off (0)
            imgArray(:,:,j) = atom2Image(:,:,1); %state 2 images in atom2Image
        end
    end
    
    
    
    if strcmp(curr,varstring) || strcmp(curr,varstring2)
        k=k+1;
        %logholdtime(ceil(i/(length(C)/TotalImages))) = str2num(next);
        logholdtime(k) = str2num(next);
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

if(sidecam)
    ftsImage = [];
    disp('Processing FTS Images...');
    
    for i=1:length(fileLocList(:))
        
        [ftsImage] = PullFTS(fileLocList{i},raw);
        imgArray(:,:,i) = ftsImage;
        
    end
    
end
    
%-------------------------------------------------------------------------%
%Average data points and calculate Centre-of-Mass movements:

%Region of Interests, use same for each image.
%Pulled from logfile:
ROIy = ROIbottom:ROItop;
ROIx = ROIleft:ROIright;

if(sidecam)
    %ROIy = 276:323; %Top:Bot
    %ROIx = 400:515;
    ROIy = 550:600; %Top:Bot
    ROIx = 520:680;
end

if(topcam)
    %ROIy = 38:88; %Top:Bot
    %ROIx = 55:110;
    ROIy = 50:100; %Top:Bot
    ROIx = 70:145;
end
center = [(ROIx(length(ROIx))-ROIx(1))/2 + ROIx(1),(ROIy(length(ROIy))-ROIy(1))/2 + ROIy(1)]; %[x y]


imgYs = []; %Array of the Y profiles (heights) to find COM movement on.
imgXs = [];
COMy = []; COMx = []; %Array of Image centre of masses of the y profiles.
avCOMy = []; avCOMx = []; %Shots per time average of centre of masses.
%image(imgArray(ROIx,ROIy,1))

xl=ROIx(1):(length(ROIx)+ROIx(1));
yl=ROIy(1):(length(ROIy)+ROIy(1));

%sum images along x direction:
for i=1:TotalImages
    for j=1:length(ROIy)
            imgYs(j,i) = sum(imgArray(ROIy(j),xl,i)); %Sum each horizontal 1D array
    end
end

%sum images along y direction:
for i=1:TotalImages
    for j=1:length(ROIx)
            imgXs(j,i) = sum(imgArray(yl,ROIx(j),i)); %Sum each vertical 1D array
    end
end

%find COM in y direction:

yvector = 1:length(ROIy);

for i=1:TotalImages
        COMy(i) = sum(imgYs(:,i).*yvector(:)) / sum(imgYs(:,i)); %Sum of Intensity*pixel location / sum of intensity
end

%find COM in x direction:

xvector = 1:length(ROIx);

for i=1:TotalImages
        COMx(i) = sum(imgXs(:,i).*xvector(:)) / sum(imgXs(:,i)); %Sum of Intensity*pixel location / sum of intensity
end

%Averaging over the three shots given that they are consecutive in the
%imgArray. Which they are not in general.

%Begin by re-arranging the COMy matrix to make contiguous with hold time:

if length(COMy)~=length(logholdtime)
    %ERROR
    error = 'error in array lengthsy'
    break
end
if length(COMx)~=length(logholdtime)
    %ERROR
    error = 'error in array lengthsx'
    break
end

COMyOrd = [];
COMxOrd = [];

%The following for loop adds another column to the COMy array which
%contains the hold time values.
for i=1:length(COMy)
    for j=1:2
        if j==1
            COMyOrd(i,j) = COMy(i);
            COMxOrd(i,j) = COMx(i);
        else
            COMyOrd(i,j) = logholdtime(i);
            COMxOrd(i,j) = logholdtime(i);
        end
    end
end

COMyOrd = sortrows(COMyOrd,2);
COMxOrd = sortrows(COMxOrd,2);%sort the rows by holdtime
    
j=1;
prevTimestep = COMyOrd(1,2);
numThisTime = 0;
sumTs = [];
avCOMy = []; avCOMx = [];
stdCOMy = []; stdCOMx = [];

%Average over the x data points per time step where x can vary from
%timestep to timestep.
for i=1:TotalImages+1
    
    if i<=TotalImages
        currTimestep = COMyOrd(i,2);
    end
    
    if currTimestep == prevTimestep && i<=TotalImages
        numThisTime = numThisTime + 1;
    else
        %new time step%
        avCOMy(j,1) = mean(sumTs(:));
        avCOMy(j,2) = prevTimestep;
        stdCOMy(j,1) = std(sumTs(:));
        stdCOMy(j,2) = prevTimestep;
        
        j=j+1;
        sumTs = [];
        numThisTime = 1;
    end

    if i<=TotalImages
        sumTs(numThisTime) = COMyOrd(i,1);
    end   

    prevTimestep = currTimestep;

end

%NOW FOR X (Horiz. direction)%%%%%%%%%%%%%%
sumTs = [];
j = 1;
numThisTime = 0;
prevTimestep = COMxOrd(1,2);
for i=1:TotalImages+1
    
    if i<=TotalImages
        currTimestep = COMxOrd(i,2);
    end
    
    if currTimestep == prevTimestep && i<=TotalImages
        numThisTime = numThisTime + 1;
    else
        avCOMx(j,1) = mean(sumTs(:));
        avCOMx(j,2) = prevTimestep;
        stdCOMx(j,1) = std(sumTs(:));
        stdCOMx(j,2) = prevTimestep;
        
        j=j+1;
        sumTs = [];
        numThisTime = 1;
    end

    if i<=TotalImages
        sumTs(numThisTime) = COMxOrd(i,1);
    end   

    prevTimestep = currTimestep;

end

avAllCOM = mean(COMyOrd(:,1)); %Centre of centre of masses
avAllCOMx = mean(COMxOrd(:,1)); %Centre of centre of masses x

%-------------------------------------------------------------------------%
%Fit periodic function to data:

%avCOMyRealtime = []; %Lets fit to real times
avAllCOMlb = avAllCOM/2;
avAllCOMub = avAllCOM*2;

avAllCOMlbx = avAllCOMx/2;
avAllCOMubx = avAllCOMx*2;

%fg = @(p,x)(p(1) + p(2)*sin(2*pi*p(3)*10^(-6).*x + p(4))) + (x*(-p(5))); %function to fit with
%p0 = [avAllCOM 1 freqguess 20 4/850];
%lb = [avAllCOMlb 0.1 5000 -100 1/500];
%ub = [avAllCOMub 5 10000 100 5];

%fg = @(p,x)(p(1) + p(2)*sin(2*pi*p(3)*10^(-6).*x + p(4))); %function to fit with
%p0 = [avAllCOM 1 freqguess 20];
%lb = [avAllCOMlb 0.1 4000 -100];
%ub = [avAllCOMub 5 10000 100];

%fg = @(p,x)(p(1) + p(2).*exp(-x./p(3)).*sin(2*pi*p(4)*10^(-6).*x + p(5))); %function to fit with damped with exp
%p0 = [avAllCOM 2 2500 freqguess 20];
%lb = [avAllCOMlb 0.2 1000 2000 -100];
%ub = [avAllCOMub 5 40000 4000 100];

fg = @(p,x)(p(1) + p(2).*exp(-x./p(3)).*sin(2*pi*p(4)*10^(-3).*x + p(5)) + (p(6)).*x); %function to fit with damped with exp
p0 = [avAllCOM 5 2000 240 20 1/200000];
lb = [avAllCOMlb 0.1 100 200 -100 1/200000];
ub = [avAllCOMub 20 2500 400 100 1/1000];

%xs = 0:timestep:(length(avCOMy)-1)*timestep;
xsY = avCOMy(:,2);
maxTimestep = max(COMyOrd(:,2));

%[coefs,r,J,cov,mse] = nlinfit(xs,avCOMy,fg,p0,'Weights',stdCOMy); %With
%weights
%[coefs,r,J,cov,mse] = nlinfit(xs,avCOMy,fg,p0);

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000);
%coefs = lsqcurvefit(fg,p0,xs,avCOMy(:,1),lb,ub,curvefitoptions); %LSQ curve fit%
[coefsY,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xsY,avCOMy(:,1),lb,ub,curvefitoptions);

fg = @(p,x)(p(1) + p(2).*exp(-x./p(3)).*sin(2*pi*p(4)*10^(-3).*x + p(5)) + (p(6)).*x); %function to fit with damped with exp
p0 = [avAllCOMx 5 2000 freqguess 20 1/200000];
lb = [avAllCOMlbx 0.1 100 10 -100 1/200000];
ub = [avAllCOMubx 20 2500 50 100 1/1000];

%xs = 0:timestep:(length(avCOMy)-1)*timestep;
xsX = avCOMx(:,2);
maxTimestep = max(COMyOrd(:,2));

%[coefs,r,J,cov,mse] = nlinfit(xs,avCOMy,fg,p0,'Weights',stdCOMy); %With
%weights
%[coefs,r,J,cov,mse] = nlinfit(xs,avCOMy,fg,p0);

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000);
%coefs = lsqcurvefit(fg,p0,xs,avCOMy(:,1),lb,ub,curvefitoptions); %LSQ curve fit%
[coefsX,resnorm,rx,exitflag,output,lambda,Jx] = lsqcurvefit(fg,p0,xsX,avCOMx(:,1),lb,ub,curvefitoptions);


coefsY
ciY = nlparci(coefsY,r,'jacobian',J)  
coefsX
ciX = nlparci(coefsX,rx,'jacobian',Jx) 
%[ypred,delta] = nlpredci(fg,xser,coefs,r,'Jacobian',J)

largexs = 0.1:0.1:maxTimestep;

%[ypred,delta] = nlpredci(fg,largexs,coefs,r,'Jacobian',J,...
%                         'MSE',mse,'SimOpt','on');
                      
% lower = ypred - delta;
% upper = ypred + delta;

%Plotting data points and fit:
%plot(largexs,fg(coefs,largexs));
%hold on;
figure(10);
errorbarxy(xsY,avCOMy(:,1),zeros(20),stdCOMy(:,1),zeros(20),stdCOMy(:,1),'.');
hold on;
plot(largexs,fg(coefsY,largexs),'k','LineWidth',0.5)
%plot(largexs,[lower;upper],'g--','LineWidth',1)
hold off;
title(['COM Vertical (y) Direction. f = ' num2str(coefsY(4))  ' (' num2str(ciY(4,1),4) ',' num2str(ciY(4,2),4) ')' ]);

figure(11);
errorbarxy(xsX,avCOMx(:,1),zeros(20),stdCOMx(:,1),zeros(20),stdCOMx(:,1),'.');
hold on;
plot(largexs,fg(coefsX,largexs),'k','LineWidth',0.5)
%plot(largexs,[lower;upper],'g--','LineWidth',1)
hold off;
title(['COM Horizon (x) Direction. f = ' num2str(coefsX(4)) ' (' num2str(ciX(4,1),4) ',' num2str(ciX(4,2),4) ')' ]);

%suggested new regions:%
SuggestedTopRegion = floor((ROIy(1) + avAllCOM) + 10)
SuggestedBottomRegion = floor((ROIy(1) + avAllCOM) - 10)


%-------------------------------------------------------------------------%
%Fancy images plot%

collated = [];

while j*length(ROIx) <= 2000
    j = j+1;
end

intWidths = j;
k=0;
total = 1;
i=1;
while total <= TotalImages
    
    if i==1
        collated(1+k*length(ROIy):(k*(length(ROIy))+length(ROIy)),...
        i:(length(ROIx)),1)...
        = imgArray(ROIy,ROIx,total);
        i = i+1;
        total = total +1;
    else
        collated(k*length(ROIy)+mod(1,k):(k*length(ROIy)+length(ROIy)-1+mod(1,k)),...
        (i-1)*length(ROIx):(i*length(ROIx)-1),1)...
        = imgArray(ROIy,ROIx,total);
        
        i = i+1;
        total = total+1;
    
        if mod(i,ceil(TotalImages/intWidths)+1) == 0
            k = k+1;         
            i = 1;
        end
    end
    
   
end

figure(2);
image(collated);

hold off;

%-------------------------------------------------------------------------%
%Fourier Transform of the data using FFT%

fgx = @(p,x)(x*-p(1)); %linear function

timestep = avCOMy(1,2); %shortest timestep of the data (that isn't zero)
if timestep == 0
    timestep = avCOMy(2,2);
end

Fs = 1/(timestep*10^(-3)); %Sampling times per second (ie once every 50us)
%Fs = 1/50;
L = max(COMyOrd(:,2)); %length of signal (maximum of the timesteps, assumes beginning at 0s

t = (0:L)*(1/Fs);
%y = 0.7*sin(2*pi*5000*t) + sin(2*pi*1200*t); %Test function
y = avCOMy(:,1);
%y = y - (fgx(coefs(4),xs) + abs(mean(fgx(coefs(4),xs)))); %For taking the
%linear out
y = y - avAllCOM;
%y = fg(coefs,t) - avAllCOM;


NFFT = 2^nextpow2(L);
%NFFT = 10000;
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(3);
plot(f,2*abs(Y(1:NFFT/2+1)));

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
    
    for j=ROIy
        for k=ROIx
    sumImages(j,k) = sum(imgArray(j,k,:));
        end
    end
    
figure(4);
image(sumImages(ROIy,ROIx),'CDataMapping','scaled');
hold off;





