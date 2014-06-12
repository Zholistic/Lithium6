%2D Trap Side Imaging Frequency Calculation Script
%Tyson Peppler jan-2014
%-------------------------------------------------------------------------%
%Load Images: This section is filename format dependent.

directory = 'C:\Data\140214_2D_Trap_Focus_Alignment_Z_4_47mm_1_75W_30us\';
datestring = '140214';
varstring = 'holdtime';
freqguess = 3000; %2pi*x Hz guess

%shotsPerTime = 4; %X points per time step% %Ver 2.0 Pull this out of
%logfile%
%times = 0:30:510; %0s to 300us in 20us steps%
%timestep = 50;

warning('OFF');
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
    
    if strcmp(curr,'ROItop')
        ROItop = next;
    elseif strcmp(curr,'ROIbottom')
        ROIbottom = next;
    elseif strcmp(curr,'ROIleft')
        ROIleft = next;
    elseif strcmp(curr,'ROIright')
        ROIright = next;
    end
                
    
end

imgArray = [];
prev = 0;
j=0;

%for loop to get information out of the logfile and load it into arrays%
for i=1:(length(C)-1)
    
    if i>1
        prev = C{i-1};
    end
    curr = C{i};
    next = C{i+1};
    %ROImax = C{i+20};
    ROImax = '1';
    
    
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
    
    if strcmp(curr,varstring)
        logholdtime(ceil(i/(length(C)/TotalImages))) = str2num(next);
    end

end


%The following for loop loads the images out of the fts files ending with
%_001.fts into 'imgArray' which is a 3d array. If images of different sizes
%use a cell data array.

%for i=1:TotalImages
    
    %filename = '140127_SIDECAM_2D_Trap_Frequency_Optical_Method';
    %filename = strcat(filename, num2str(times(floor(((i-1)+shotsPerTime)/shotsPerTime)))); %ie times(1) will be 0 for the first 5 iterations%    
    
    %midString = 'Point  ';
    %filename = [filename midString num2str(mod(i,shotsPerTime)+1)];
    
    %if i<10
    %  picno = ['_00' num2str(i) '_00' num2str(i)];
    %else
    %  picno = ['_0' num2str(i) '_0' num2str(i)];  
    %end
    %picno
    
    %fullname = [directory filename picno '.fts'];
    %fullname = '140110_2DTrap_5W_TrapFrequency_2msTOF_250us_Release_Recapture_Constant_gradient_check_shot_to_shot_COM_10us_ImagePulse_GradStayOff_HoldTime_0Point  1_001.fts'
    
    
    %ftsImg = imread(fullname); %ImgArray will be the array of all images loaded with load_fts%
    %imgArray(:,:,i) = ftsImg; %image(imgArray(1360x,1024y,TotalImages)) access format with sizes
    
    %TODO: Try difference between atom and beam images instead of the
    %_001.fts image.
    
    
%end
    

%-------------------------------------------------------------------------%
%Average data points and calculate Centre-of-Mass movements:

%Region of Interests, use same for each image.
ROIy = str2num(ROIbottom):str2num(ROItop);
ROIx = str2num(ROIleft):str2num(ROIright);

ROIy = 276:323; %Top:Bot
ROIx = 400:515;

center = [(ROIx(length(ROIx))-ROIx(1))/2 + ROIx(1),(ROIy(length(ROIy))-ROIy(1))/2 + ROIy(1)]; %[x y]

imgYs = []; %Array of the Y profiles (heights) to find COM movement on.
COMy = []; %Array of Image centre of masses of the y profiles.
avCOMy = []; %Shots per time average of centre of masses.

%image(imgArray(ROIx,ROIy,1))

xl=ROIx(1):(length(ROIx)+ROIx(1));


%sum images along x direction:
for i=1:TotalImages
    for j=1:length(ROIy)
            imgYs(j,i) = sum(imgArray(ROIy(j),xl,i)); %Sum each horizontal 1D array
    end
end

%find COM in y direction:

yvector = 1:length(ROIy);

for i=1:TotalImages
        COMy(i) = sum(imgYs(:,i).*yvector(:)) / sum(imgYs(:,i)); %Sum of Intensity*pixel location / sum of intensity
end

%Averaging over the three shots given that they are consecutive in the
%imgArray. Which they are not in general.

%Begin by re-arranging the COMy matrix to make contiguous with hold time:

if length(COMy)~=length(logholdtime)
    %ERROR
    error = 'error in array lengths'
end

COMyOrd = [];

%The following for loop adds another column to the COMy array which
%contains the hold time values.
for i=1:length(COMy)
    for j=1:2
        if j==1
            COMyOrd(i,j) = COMy(i);
        else
            COMyOrd(i,j) = logholdtime(i);
        end
    end
end

COMyOrd = sortrows(COMyOrd,2); %sort the rows by holdtime
    
j=1;
prevTimestep = COMyOrd(1,2);
numThisTime = 0;
sumTs = [];
avCOMy = [];
stdCOMy = [];

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

%errorbarxy(1:16,avCOMy,zeros(16),stdCOMy,zeros(16),stdCOMy,'.')

avAllCOM = mean(COMyOrd(:,1)); %Centre of centre of masses

%-------------------------------------------------------------------------%
%Fit periodic function to data:

%avCOMyRealtime = []; %Lets fit to real times
avAllCOMlb = avAllCOM/2;
avAllCOMub = avAllCOM*2;

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

fg = @(p,x)(p(1) + p(2).*exp(-x./p(3)).*sin(2*pi*p(4)*10^(-6).*x + p(5)) + (p(6)).*x); %function to fit with damped with exp
p0 = [avAllCOM 2 2500 freqguess 20 1/200000];
lb = [avAllCOMlb 1.5 1963 2000 -100 1/200000];
ub = [avAllCOMub 5 10468 8000 100 1/1000];

%xs = 0:timestep:(length(avCOMy)-1)*timestep;
xs = avCOMy(:,2);
maxTimestep = max(COMyOrd(:,2));

%[coefs,r,J,cov,mse] = nlinfit(xs,avCOMy,fg,p0,'Weights',stdCOMy); %With
%weights
%[coefs,r,J,cov,mse] = nlinfit(xs,avCOMy,fg,p0);

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000);
%coefs = lsqcurvefit(fg,p0,xs,avCOMy(:,1),lb,ub,curvefitoptions); %LSQ curve fit%
[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,avCOMy(:,1),lb,ub,curvefitoptions);

coefs
%[ypred,delta] = nlpredci(fg,xser,coefs,r,'Jacobian',J)

largexs = 0.1:0.1:maxTimestep;

%[ypred,delta] = nlpredci(fg,largexs,coefs,r,'Jacobian',J,...
%                         'MSE',mse,'SimOpt','on');
  

ci = nlparci(coefs,r,'jacobian',J)               
                     
% lower = ypred - delta;
% upper = ypred + delta;

%Plotting data points and fit:
%plot(largexs,fg(coefs,largexs));
%hold on;
figure(1);
errorbarxy(xs,avCOMy(:,1),zeros(20),stdCOMy(:,1),zeros(20),stdCOMy(:,1),'.');
hold on;
plot(largexs,fg(coefs,largexs),'k','LineWidth',0.5)
%plot(largexs,[lower;upper],'g--','LineWidth',1)
hold off;

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

Fs = 1/(timestep*10^(-6)); %Sampling times per second (ie once every 50us)
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
image(sumImages(ROIy,ROIx));
hold off;





