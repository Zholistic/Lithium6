function [ ] = EOSGenerate(inputDensity, inputRadius, calcRegion, omegaR)

massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
a0 = 5.29e-11; %o.0

pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
kpixelLength = (13e-6*(83/400)); %kristian pixel length

%radproftoFityReal = radProfilesAvg(1,:,1)./(kpixelLength^2); %convert to density/m^-2
radproftoFityReal = inputDensity; %density already in m^-2
radiusVector = inputRadius; %x vector 
radiusVectorMeters = radiusVector.*kpixelLength; %convert x vector to m
potential = 0.5 .* massL6 .* omegaR^2 .* radiusVectorMeters.^2;
potential_nk = (potential./kB).*(10^9);
%plot(potential_nk,radProfilesAvg(1,:,i));

corrFactor = 1;
radProfilesDensityAdjusted = radproftoFityReal.*corrFactor; %2 for spin states and 1.2 for correction factor

a0 = 5.29e-11; %o.0
a2d = 301946 * a0; %972G, 5.789*2pi hz
eb = (hbar^2)/(a2d^2 * massL6);


%------------kappa and p calculations:
%calcRegion = 2:650;
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

close all;
%averaged cloud kappa and p's
figure(8); 
subplot(3,2,1); plot(radiusVectorMeters,radProfilesDensityAdjusted,'.');
title('Density Adjusted vs radius (meters)');
subplot(3,2,5); plot(pTildeAvg(3:end),kappaTildeAvg(3:end),'.');
grid on;
title('kappaTilde vs pTilde'); axis([0 14 0 3]);
subplot(3,2,3); plot(potential(calcRegion),pTildeAvg,'.');
grid off;
title('pTilde vs potential (J)');
subplot(3,2,4); plot(potential(calcRegion(3:end)),kappaTildeAvg(3:end),'.');

title('kappaTilde vs potential (J)');
subplot(3,2,2); plot(potential(calcRegion),radProfilesDensityAdjusted(calcRegion),'.');
title('density vs potential (J)');

%catchArray = [];

    %catchArray(1,:) =[pTildeAvg(3:end) kappaTildeAvg(3:end)];


%outputArray = catchArray;

end