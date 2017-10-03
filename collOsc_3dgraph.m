
close all;
breathingData = csvimport('C:\Users\tpeppler\Dropbox\PhD\2D_2017\ColOsc Data Review\collOscDatacellv10.csv');

breathingFreqsCell = breathingData(2:end,6);
breathingFreqErrorsCell = breathingData(2:end,7);
breathingAtomNumsCell = breathingData(2:end,8);
breathingMagfieldsCell = breathingData(2:end,3);

breathingFreqs = []; breathingFreqErrors = []; breathingAtomNums = []; breathingMagfields = [];
j=1;
for i=1:length(breathingFreqsCell)
    
    if(i <= 39 && i >= 29)
    %if(breathingFreqErrorsCell{i} < 0.7 && breathingAtomNumsCell{i} < 10E3)
    breathingFreqs(j) = breathingFreqsCell{i};
    breathingFreqErrors(j) = breathingFreqErrorsCell{i};
    breathingAtomNums(j) = breathingAtomNumsCell{i};
    breathingMagfields(j) = breathingMagfieldsCell{i};
    j = j+1;
    end
    
end

figure(10);
plot(breathingAtomNums,breathingFreqs,'.'); grid on; axis([0 12E4 35 60]);
%errorbar(breathingAtomNums,breathingFreqs,breathingFreqErrors,breathingFreqErrors,'.'); grid on; axis([0 12E4 35 60]);

figure(11);
%plot(breathingMagfields,breathingFreqs,'.'); grid on; axis([650 900 35 60]);
errorbar(breathingMagfields,breathingFreqs,breathingFreqErrors,breathingFreqErrors,'.'); grid on; axis([650 900 35 60]);

x = []; y = []; z = []; xlin = []; ylin = []; Z = []; X = []; Y = [];

x = breathingMagfields;
y = breathingAtomNums;
z = breathingFreqs;

xlin = linspace(min(x),max(x),50);
ylin = linspace(min(y),max(y),50);
        

[X,Y] = meshgrid(xlin,ylin);

Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear

figure(333);
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
%surf(X,Y,Z);
scatter3(x,y,z);
axis([650 900 0 12E4 35 60]); 
grid on;
xlabel('Mag Field') % x-axis label
ylabel('Atom Num') % y-axis label
zlabel('Breathing Freq / Dipole Freq') % y-axis label
view(3);



%------------Dipole section


dipoleData = csvimport('C:\Users\tpeppler\Dropbox\PhD\2D_2017\ColOsc Data Review\collOscDipoleDatacellv18.csv');

recent = 0;
dipoleFreqsCell = dipoleData(2+recent:end,6);
dipoleFreqErrorsCell = dipoleData(2+recent:end,7);
dipoleAtomNumsCell = dipoleData(2+recent:end,8);
dipoleMagfieldsCell = dipoleData(2+recent:end,3);
dipoleFreqsxCell = dipoleData(2+recent:end,9);
dipoleFreqsxErrorsCell = dipoleData(2+recent:end,10);
dipoleFreqsyCell = dipoleData(2+recent:end,11);
dipoleFreqsyErrorsCell = dipoleData(2+recent:end,12);
dipoleFreqGeoMeanCell = dipoleData(2+recent:end,13);
dipoleFreqPCACell = dipoleData(2+recent:end,14);
dipoleFreqPCAErrorsCell = dipoleData(2+recent:end,15);
dipoleAspectRatioCell = dipoleData(2+recent:end,16);

dipoleFreqs = []; dipoleFreqErrors = []; dipoleAtomNums = []; dipoleMagfields = []; dipoleFreqsx = []; dipoleFreqsyError = [];
dipoleFreqsGeoMean = []; dipoleFreqsy = []; dipoleFreqsPCA = []; dipoleFreqsPCAError = []; dipoleAspectRatio = [];
dipoleFreqXdivY = []; dipoleFreqXminusY = [];
j=1;
for i=1:length(dipoleFreqsCell)
    
    %constrain the points to ones where the x&y error on the frequency is
    %less than 0.1
    if(dipoleFreqsxErrorsCell{i} < 0.3 && dipoleFreqsyErrorsCell{i} < 0.3 && dipoleAtomNumsCell{i} < 60E3)
    %if(dipoleFreqPCAErrorsCell{i} < 0.2 && dipoleAtomNumsCell{i} < 70E3)
    dipoleFreqs(j) = dipoleFreqsCell{i};
    dipoleFreqErrors(j) = dipoleFreqErrorsCell{i};
    dipoleAtomNums(j) = dipoleAtomNumsCell{i};
    dipoleMagfields(j) = dipoleMagfieldsCell{i};
    dipoleFreqsx(j) = dipoleFreqsxCell{i};
    dipoleFreqsxError(j) = dipoleFreqsxErrorsCell{i};
    dipoleFreqsy(j) = dipoleFreqsyCell{i};
    dipoleFreqsyError(j) = dipoleFreqsyErrorsCell{i};
    dipoleFreqsGeoMean(j) = dipoleFreqGeoMeanCell{i};
    dipoleFreqsPCA(j) = dipoleFreqPCACell{i};
    dipoleFreqsPCAError(j) = dipoleFreqPCAErrorsCell{i};
    dipoleAspectRatio(j) = dipoleAspectRatioCell{i};
    dipoleFreqXminusY(j) = dipoleFreqsx(j) - dipoleFreqsy(j);
    dipoleFreqXdivY(j) = dipoleFreqsx(j)./dipoleFreqsy(j);
    j = j+1;
    end
   
end



%average the dipole frequency at each magnetic field:
sortedDipoleMagfields = []; sortedDipoleFreqsx = []; sortedDipoleFreqsy = []; sortedDipoleFreqsGeoMean = []; extantMagfields = [];
sortedDipoleFreqsPCA = [];

[sortedDipoleMagfields,indexs] = sort(dipoleMagfields);
sortedDipoleFreqsx = dipoleFreqsx(indexs);
sortedDipoleFreqsy = dipoleFreqsy(indexs);
sortedDipoleFreqsGeoMean = dipoleFreqsGeoMean(indexs);
sortedDipoleFreqsPCA = dipoleFreqsPCA(indexs);

extantMagfields = unique(sortedDipoleMagfields);

freqx_avg_at_fields = []; freqy_avg_at_fields = []; freqGeoMean_avg_at_fields = []; freqPCA_avg_at_fields = [];
for i=1:length(extantMagfields)
    currField = extantMagfields(i);
    
    freqx_avg_at_fields(i) = mean(abs(sortedDipoleFreqsx(sortedDipoleMagfields == currField)));
    freqy_avg_at_fields(i) = mean(abs(sortedDipoleFreqsy(sortedDipoleMagfields == currField)));
    freqGeoMean_avg_at_fields(i) = mean(abs(sortedDipoleFreqsGeoMean(sortedDipoleMagfields == currField)));
    freqPCA_avg_at_fields(i) = mean(abs(sortedDipoleFreqsPCA(sortedDipoleMagfields == currField)));
    
end 

%Plots:

figure(21);
plot(dipoleMagfields,dipoleFreqsGeoMean,'.'); grid on; axis([650 900 20 30]); hold on;
plot(extantMagfields,freqGeoMean_avg_at_fields,'x');
xss = 0:1:900; plot(xss,1.068*sqrt(xss) + (-5.639),'-.');
%v14, 1.19 bounds(0.975, 1.4) * sqrt(x) + (-8.921 bounds(-14.79, -3.049))
%v14, 0.02139*x
%plot(xss,0.975*sqrt(xss) - 8.921,'-');
%plot(xss,1.4*sqrt(xss) - 8.921,'-');
%v15 below:
% f(x) = a*x^(1/2)+b
%Coefficients (with 95% confidence bounds):
%       a =       1.068  (0.9463, 1.19)
%       b =      -6.639  (-8.963, -2.314)\
%v16, below:
%General model:
%     f(x) = a*x^(1/2)+b
%Coefficients (with 95% confidence bounds):
%       a =       1.068  (0.9463, 1.19)
%       b =      -5.639  (-8.963, -2.314)
       
figure(22);
plot(dipoleMagfields,dipoleFreqsx,'.'); grid on; axis([650 900 10 35]); 
hold on;
plot(extantMagfields,freqx_avg_at_fields,'xb');
plot(dipoleMagfields,dipoleFreqsy,'.r'); 
plot(extantMagfields,freqy_avg_at_fields,'xr');

figure(23);
plot(dipoleMagfields,dipoleFreqsPCA,'.'); grid on; axis([650 900 10 35]); hold on;
plot(extantMagfields,freqPCA_avg_at_fields,'x');

figure(24);
plot(dipoleMagfields,dipoleAspectRatio,'.'); grid on; axis([650 900 0.9 1.3]); hold on;

figure(25);
plot(dipoleMagfields,dipoleFreqXminusY,'.'); grid on; axis([650 900 -1.5 1.5]); hold on;

figure(26);
plot(dipoleAspectRatio,dipoleFreqXminusY,'.'); grid on; axis([0.95 1.05 -1.5 1.5]); hold on;

figure(27);
plot(dipoleAspectRatio,dipoleFreqXdivY,'.'); grid on; axis([0.95 1.05 0.95 1.05]); hold on;



x2 = []; y2 = []; z2 = []; xlin2 = []; ylin2 = []; Z2 = []; X2 = []; Y2 = [];

x2 = dipoleMagfields;
y2 = dipoleAtomNums;
z2 = dipoleFreqsPCA;

xlin2 = linspace(min(x2),max(x2),50);
ylin2 = linspace(min(y2),max(y2),50);
        

[X2,Y2] = meshgrid(xlin2,ylin2);

Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear

figure(366);
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
%surf(X,Y,Z);
scatter3(x2,y2,z2);
axis([650 900 0 12E4 20 30]); 
grid on;
xlabel('Mag Field') % x-axis label
ylabel('Atom Num') % y-axis label
zlabel('Dipole Freq') % y-axis label
view(3);


%function results:
%v14, 1.19 bounds(0.975, 1.4) * sqrt(x) + (-8.921 bounds(-14.79, -3.049))
%v15, 1.068*sqrt(extantMagfields(i)) - 5.639;
    %   a =      0.8765  (0.5511, 1.202)
    %   b =     -0.5441  (-9.212, 8.124)

freq_from_fit = [];
for i=1:length(extantMagfields)
    freq_from_fit(i) = 1.068*sqrt(extantMagfields(i)) -5.639;
end

dipoleFreqForBreathingFields = [];
for i=1:length(breathingMagfields)
    dipoleFreqForBreathingFields(i) = 1.068*sqrt(breathingMagfields(i)) -5.639;
end

figure(101);
plot(breathingAtomNums,breathingFreqs./dipoleFreqForBreathingFields,'.'); grid on; axis([0 12E4 0 3]);
%errorbar(breathingAtomNums,breathingFreqs,breathingFreqErrors,breathingFreqErrors,'.'); grid on; axis([0 12E4 35 60]);

figure(102);
%plot(breathingMagfields,breathingFreqs./dipoleFreqForBreathingFields,'.'); 
percentErrorBreathing = breathingFreqErrors./breathingFreqs;
errorbar(breathingMagfields,breathingFreqs./dipoleFreqForBreathingFields,(breathingFreqs./dipoleFreqForBreathingFields).*percentErrorBreathing,(breathingFreqs./dipoleFreqForBreathingFields).*percentErrorBreathing,'.'); grid on; axis([650 900 35 60]);
grid on; axis([650 900 1.8 2.3]);

x = []; y = []; z = []; xlin = []; ylin = []; Z = []; X = []; Y = [];

x = breathingMagfields;
y = breathingAtomNums;
z = breathingFreqs./dipoleFreqForBreathingFields;

xlin = linspace(min(x),max(x),50);
ylin = linspace(min(y),max(y),50);
        

[X,Y] = meshgrid(xlin,ylin);

Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear

figure(333);
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
%surf(X,Y,Z);
scatter3(x,y,z);
axis([650 900 0 12E4 0 3]); 
grid on;
xlabel('Mag Field') % x-axis label
ylabel('Atom Num') % y-axis label
zlabel('Breathing Freq') % y-axis label
view(3);

%breathingAtomNums(find(breathingAtomNums < 10000))
binnedResult = binMe(breathingFreqs(find(breathingMagfields < 705))./dipoleFreqForBreathingFields(find(breathingMagfields < 705)),breathingAtomNums(find(breathingMagfields < 705)),100);
figure(103);
plot(binnedResult(2,:),binnedResult(1,:),'.'); grid on; axis([0 60000 1 3]);

%------------------------
%trajectory plots:
if(0)
close all;
for i=1:length(collOscDipoleData)
    ys = collOscDipoleData{i}.centersy;
    xs = collOscDipoleData{i}.centersx;
    times = collOscDipoleData{i}.centersxs;
    
    figure(3000+i); plot(xs,ys); hold on; plot(xs,ys,'.');
end
end





