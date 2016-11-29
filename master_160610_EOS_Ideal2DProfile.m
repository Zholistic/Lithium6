

%omegaR = 2*pi*25; %radial trap frequency
omegaRVector = 2.*pi.*[24.88224858 25.09885405 25.66767827 26.3891196];
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
mu = kB*191e-9; %Chemical Potential 
T = 25e-9; %Temperature 25nK
%T = 0;
kpixelLength = (13e-6*(83/400)); %kristian pixel length

xValues = []; yValues = []; zValues = [];

xValues = (-100:100).*kpixelLength;
yValues = (-100:100).*kpixelLength;

lambda = sqrt((2*pi*hbar^2)/(massL6*kB*T));

%zValues = 1/(lambda^2) .* (log(1 + exp((mu - (1/2)*massL6*(omegaR^2)*(xValues.^2) - (1/2)*(omegaR^2)*(yValues.^2))/(kB * T))));
%zValues = (log(1 + exp((mu - (1/2)*massL6*(omegaR^2)*(xValues.^2) - (1/2)*(omegaR^2)*(yValues.^2))/(kB * T))));

for i=1:length(xValues)
    for j=1:length(yValues)
        zValues(i,j) = (1) .* 1/(lambda^2) .* (log(1 + exp((mu - (1/2).*massL6.*(omegaRVector(4).^2).*(xValues(i).^2) - (1/2).*massL6.*(omegaRVector(4).^2)*(yValues(j).^2))/(kB * T))));
    end
end

close all;
h = surf(zValues);
 set(h,'edgecolor','none');
 
 imageArrayAvgs = [];
 imageArrayAvgs(:,:,1) = zValues;
 
%%%%%Atom numbers:
pixelCounts = [];
for i=1:length(imageArrayAvgs(1,1,:))
    pixelCounts(i) = sum(sum(imageArrayAvgs(:,:,i)));
end

atomNumber = pixelCounts*(kpixelLength.^2);
 
%Non-radially averaged profile:
nonRadProfile = zValues(101,101:end);

%Radially averaged average profiles:
radProfilesAvgReorder = []; radProfilesTAvg = []; centerAvgs = [];
disp('Radially averaging...');
for i=1:length(imageArrayAvgs(1,1,:))
    [radProfilesTAvg(:,:,i),centerAvgs(:,i)] = radAverageReorder(imageArrayAvgs(:,:,i));
    radProfilesAvgReorder(:,:,i) = radProfilesTAvg(:,1:end-5,i);
end

radProfilesAvg = [];
radProfilesAvg(2,:,1) = meanNelements(radProfilesAvgReorder(2,:,1),20);
radProfilesAvg(1,:,1) = meanNelements(radProfilesAvgReorder(1,:,1),20);

binnedOutput = binMe(radProfilesAvg(1,:,1),radProfilesAvg(2,:,1),140);

%radProfilesAvg = [];
%radProfilesAvg(2,:,1) = binnedOutput(2,:);
%radProfilesAvg(1,:,1) = binnedOutput(1,:);

%Shift wings of radial profiles to zero:
if(0)
shiftByAvgs = [];
for i=1:length(radProfilesAvg(1,1,:))
    shiftByAvgs(i) = mean(radProfilesAvg(1,72:82,i));
    radProfilesAvg(1,:,i) = radProfilesAvg(1,:,i) - shiftByAvgs(i);
end

%Shift images by found factor:
imageArrayAvgsZerod = [];
for i=1:length(imageArrayAvgs(1,1,:))
    imageArrayAvgsZerod(:,:,i) = imageArrayAvgs(:,:,i) - shiftByAvgs(i);
end
end

%Azimuthually averaged profile vs non-averaged profile:
%plot(nonRadProfile); hold on; plot(radProfilesAvg(2,:,1),radProfilesAvg(1,:,1));

if(0)
    for i=1:length(imageArrayAvgs(1,1,:))
        if(mod(i,1) == 0)
            figure(i);
            hold on; plot(radProfilesAvg(2,:,i),radProfilesAvg(1,:,i),'r'); line([0 100],[0 0]); hold off;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%----KappaTilde vs PTilde
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
magfield = 865;
magfield = 972;
omegaz = 5789 * 2 * pi;
omegar = 25 * 2 * pi; %Double check this weak trapping


a0 = 5.29e-11; %o.0
a2dVector = a0.*[47194.1 66130.5 141838 302318];
magVector = [865 880 920 972];
az = sqrt(hbar/(massL6*omegaz));

pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
kpixelLength = (13e-6*(83/400)); %kristian pixel length

%radproftoFityReal = radProfilesAvg(1,:,1)./(kpixelLength^2); %convert to density/m^-2
radproftoFityReal = radProfilesAvg(1,:,1); %density already in m^-2
radiusVector = radProfilesAvg(2,:,1); %x vector 
radiusVectorMeters = radiusVector.*kpixelLength; %convert x vector to m
potential = 0.5 .* massL6 .* omegaRVector(4)^2 .* radiusVectorMeters.^2;
potential_nk = (potential./kB).*(10^9);
%plot(potential_nk,radProfilesAvg(1,:,i));

%convert radProfile pixels to potential:
%radProfilesPotential = 0.5 .* massL6 .* omegaRVector(4)^2 .*(radProfilesAvg(2,:,1).*kpixelLength).^2;
%plot(radProfilesPotential,radProfilesAvg(1,:,1)./(pixelLength^2));
%radProfilesPotentialConv = radProfilesPotential.*10^30;
radProfilesDensityAdjusted = radproftoFityReal.*2; %2 for spin states and 1.2 for correction factor
%plot(potential_nk,radProfilesDensityAdjusted);
%plot(potential,radproftoFityReal);

a0 = 5.29e-11; %o.0
a2d = 47186.7 * a0; %From mathematica spreadsheet TODO generate text file
a2d = 301946 * a0; %972G, 5.789*2pi hz
eb = (hbar^2)/(a2d^2 * massL6);


%Scale x&y to H.O. dimensionless 
%Density SI = m^-2
radProfilesDensityConv = radProfilesDensityAdjusted./(massL6*omegaz / hbar);
%Potential SI = kg * s^-2 * m^2 // Energy units = kg * m^2 * s^-2
radProfilesPotentialConv = potential./((hbar*omegaz)/2); 
% plot(radProfilesPotentialConv,radProfilesDensityConv);


%------------kappa and p calculations:
calcRegion = 2:650;
kappaTildeAvg = []; pTildeAvg = []; pAvg = []; dydx = [];

dydx = diff([eps radProfilesDensityAdjusted(calcRegion)])./diff([eps potential(calcRegion)]);
%kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) * gradient(radProfilesDensityAdjusted(calcRegion),potential(calcRegion));
kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) .* dydx;

for j=1:length(calcRegion)
    pAvg(j) = trapz(potential(j:calcRegion(end)),radProfilesDensityAdjusted(j:calcRegion(end))); %integral from end to j point
    n = radProfilesDensityAdjusted(j);     %n, density on j'th point (varies across the potential)
    pTildeAvg(j) = ((2*massL6)/(n^2 * hbar^2 * pi)) * pAvg(j);
end
%-------------

%------------kappa and p HO units calculations:
%calcRegion = 7:10000;
kappaTildeAvgHO = []; pTildeAvgHO = []; pAvgHO = []; dydx = [];

dydx = diff([eps radProfilesDensityConv(calcRegion)])./diff([eps radProfilesPotentialConv(calcRegion)]);
%dydx convert from HO units
dydx(2:end) = dydx(2:end) .*(massL6*omegaz / hbar) .* (1/(((hbar*omegaz)/2)));
%kappaTildeAvgHO = (-1)*(pi * hbar^2 / massL6) * gradient(radProfilesDensityConv(calcRegion),radProfilesPotentialConv(calcRegion)).*(massL6*omegaz / hbar) .* (1/(((hbar*omegaz)/2)));
kappaTildeAvgHO = (-1)*(pi * hbar^2 / massL6) .* dydx;

for j=1:length(calcRegion)
    pAvgHO(j) = trapz(radProfilesPotentialConv(j:calcRegion(end)),radProfilesDensityConv(j:calcRegion(end))); %integral from end to j point
    n = radProfilesDensityConv(j).*(massL6*omegaz / hbar);    %n, density on j'th point (varies across the potential)
    pAvgHO(j) = pAvgHO(j).*((hbar*omegaz)/2).*(massL6*omegaz / hbar);
    pTildeAvgHO(j) = ((2*massL6)/(n^2 * hbar^2 * pi)) * pAvgHO(j);
end
%-------------

%Scale x&y from H.O. dimensionless 
profDensity = []; profPotential = [];
%Density SI = m^-2
profDensity = radProfilesDensityConv.*(massL6*omegaz / hbar);
%Potential SI = kg * s^-2 * m^2 // Energy units = kg * m^2 * s^-2
profPotential = potential.*((hbar*omegaz)/2); 
% plot(radProfilesPotentialConv,radProfilesDensityConv);

close all;
%averaged cloud kappa and p's
figure(8); 
subplot(3,2,1); plot(radiusVectorMeters,radProfilesDensityAdjusted,'.');
title('Density Adjusted vs radius (meters)');
subplot(3,2,5); plot(pTildeAvg(3:end),kappaTildeAvg(3:end),'.');
grid on;
title('kappaTilde vs pTilde'); axis([0 14 0 2]);
subplot(3,2,3); h1 = plot(potential(calcRegion),pTildeAvg,'.');
grid off;
title('pTilde vs potential (J)');
subplot(3,2,4); plot(potential(calcRegion(3:end)),kappaTildeAvg(3:end),'.');

title('kappaTilde vs potential (J)');
subplot(3,2,2); plot(potential(calcRegion),radProfilesDensityAdjusted(calcRegion),'.');
title('density vs potential (J)');

%HO units:
%close all;
%averaged cloud kappa and p's
%figure(88); plot(radiusVectorMeters,radProfilesDensityAdjusted);
%title('Density Adjusted vs radius (meters)');
%figure(89); plot(pTildeAvgHO,kappaTildeAvgHO,'.');
%title('kappaTilde vs pTilde'); axis([0 14 0 2]);
%figure(90); h1 = plot(potential(calcRegion),pTildeAvgHO,'.');
%title('pTilde vs potential (J)');
%figure(91); plot(potential(calcRegion),kappaTildeAvgHO,'.');
%title('kappaTilde vs potential (J)');
%figure(92); plot(radProfilesPotentialConv(calcRegion),radProfilesDensityConv(calcRegion),'.');
%title('HO Units: density vs potential');



