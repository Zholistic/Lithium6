directory = 'C:\Data\160211_2D_Trap_Freq_Measure\'; %4.8kHz
directory = 'C:\Data\160217_2D_Trap_Freq_Measure_14W\'; %5.32kHz
directory = 'C:\Data\160218_2D_TrapFreq_Measures\'; %5.8kHz
%directory = 'C:\Data\160205_VerdiG_Tight_Trapping_Frequency_Measurement\';

date = '160211';
date = '160217';
date = '160218';
%date = '160205';
camera = 'sidecam';
varstring = 'holdtime';

pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 134*10; %135*x us
%Isat = 10^6;
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

holdtimes = []; centerys = [];
holdtimes = varData(:,1);
centerys = varData(:,6);


centeryAvgs = []; holdtimeAvgs = []; centeryStds = []; 
k = 0; prev = 0; j = -1;
for i=1:length(centerys)
    curr = holdtimes(i);
    
    if(curr==prev)
        j = j+1;
        
    else
        k=k+1;
        holdtimeAvgs(k) = mean(holdtimes(i-j-1:i-1));
        centeryAvgs(k) = mean(centerys(i-j-1:i-1));
        centeryStds(k) = std(centerys(i-j-1:i-1));
        j=0;
        
    end
    prev = curr;
    
end

if(0)
plot(holdtimes,centerys,'.');
hold on;
plot(holdtimeAvgs,centeryAvgs,'.r');
errorbar(holdtimeAvgs,centeryAvgs,centeryStds/2,'.r');
end

centeryAvgsZerod = [];
for i=1:length(centeryAvgs)
centeryAvgsZerod(i) = centeryAvgs(i) - mean(centeryAvgs);
end

holdtimeAvgsSeconds = holdtimeAvgs.*10^(-6);

centeryAvgsFirstSeq = []; holdtimeAvgsSecondsFirstSeq = [];
holdtimeAvgsSecondsFirstSeq = holdtimeAvgsSeconds(1:16);
centeryAvgsFirstSeq = centeryAvgs(1:16);

centerAvgs160218 = [centeryAvgs(4:21) centeryAvgs(35:end)];
holdtimeAvgsSeconds160218 = [holdtimeAvgsSeconds(4:21) holdtimeAvgsSeconds(35:end)];
centerStds160218 = [centeryStds(4:21) centeryStds(35:end)];

sineExpDamp = @(p,x)(p(1).*exp(-p(2).*x).*sin(p(3).*x+p(4))-p(5)*x+p(6));
gcoefsSineDamp = []; gcoefsSineDampError = [];
[gcoefsSineDamp, gcoefsSineDampError] = sinExpDampLinModFit(centeryAvgs,holdtimeAvgsSeconds); 

plot(sineExpDamp(gcoefsSineDamp,0:0.000001:0.004),'r');
hold on;
errorbar(holdtimeAvgsSeconds./(10^(-6)),centeryAvgs,centeryStds/2,'.');

freq = gcoefsSineDamp(3)/(2*pi);





