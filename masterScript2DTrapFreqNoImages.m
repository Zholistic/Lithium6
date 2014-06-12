directory = 'C:\Data\140605_2DTrapAlignment_5WattInCell_Sidecam\';
date = '140605';
camera = 'sidecam';
varstring = 'holdtime';
varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %2.84 um topcam
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 135*10; %135*x us
kB = 1.38e-23; %Boltzmanns constant m^2 kg s^-2 K^-1
imgArrayFresh = [];  lowIntRealAtomImg = [];
OD = 0; %optical density from SPE process function 1=OD, 0=WithSigma
close all;


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
    
    
    
    if strcmp(curr,varstring) || strcmp(curr,varstring2)
        k=k+1;
        varData(k,1) = str2num(next);
    end
    
    %Also pull centers out:
    if strcmp(curr,'centery3W')
        j=j+1;
        varData(j,2) = str2num(next);
    elseif strcmp(curr,'centery')
        j=j+1;
        varData(j,2) = str2num(next);
    elseif strcmp(curr,'centery2')
        j=j+1;
        varData(j,2) = str2num(next);
    elseif strcmp(curr,'centery3')
        j=j+1;
        varData(j,2) = str2num(next);
    end
    
       

end


%All in varData:
if(0)
threeWattData = [];
threeWattData = varData(63:end,:);
sData = sortrows(threeWattData(:,:));
timestep = 25;
avgData = []; stdDevData = [];

j=0;
for i=1:3:length(sData)
    j=j+1;
    avgData(j,2) = mean(sData(i:i+2,2));
    avgData(j,1) = sData(i,1);
    stdDevData(j,2) = std(sData(i:i+2,2));
    stdDevData(j,1) = sData(i,1);
end
end
fiveWattData = [];
fiveWattData = varData(27:end,:);
sData = sortrows(fiveWattData(:,:));
timestep = 25;
avgData = []; stdDevData = [];

j=0;
for i=1:3:length(sData)-3
    j=j+1;
    avgData(j,2) = mean(sData(i:i+2,2));
    avgData(j,1) = sData(i,1);
    stdDevData(j,2) = std(sData(i:i+2,2));
    stdDevData(j,1) = sData(i,1);
end
    

%plot(threeWattData(:,1),threeWattData(:,2),'x');


%fg = @(p,x)(p(1) + p(2).*exp(-x./p(3)).*sin(2*pi*p(4).*x + p(5)) + (p(6)).*x); %function to fit with damped with exp
%p0 = [583 3 0.001 5500 20 890];
%lb = [570 0.1 0.00001 5000 -50 100];
%ub = [600 10 0.1 6000 50 10000];

fg = @(p,x)(p(1) + p(2).*sin(2*pi*p(3).*x + p(4)) + (p(5)).*x); %function to fit with damped with exp
p0 = [583 3 5500 20 890];
lb = [570 0.4 5000 -50 1];
ub = [600 10 8000 50 10000];

%xs = 0:timestep:(length(threeWattData)-1)*timestep;
xs = 0:10:max(avgData(:,1));

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000);
%coefs = lsqcurvefit(fg,p0,xs,avCOMy(:,1),lb,ub,curvefitoptions); %LSQ curve fit%
[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,avgData(:,1)./1e6,avgData(:,2),lb,ub,curvefitoptions);

coefs
ci = nlparci(coefs,r,'jacobian',J);  


figure(10);
errorbarxy(avgData(:,1)./1e6,avgData(:,2),zeros(20),stdDevData(:,2),zeros(20),stdDevData(:,2),'.');
hold on;
plot(xs./1e6,fg(coefs,xs./1e6),'k','LineWidth',0.5);




