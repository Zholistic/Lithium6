function [xarray, yarray ] = fitVirial(inputDensityArray, inputRadiusArray, calcRegion, omegaz, omegar, field, a2d, smoothOn, zeroOn, cutoffLeft, cutoffRight)

massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
a0 = 5.29e-11; %o.0
bfield = field;
kpixelLength = (13e-6*(83/400)); %kristian pixel length
smoothAmount = 5;
eb = (hbar^2)/(a2d^2 * massL6);


%---Zeroing and Smoothing:
inputDensityArrayS = []; shiftByThis = []; shiftByThisMin = [];
if(smoothOn)
    for i=1:length(inputDensityArray(1,1,:))
        inputDensityArrayS(1,calcRegion(1)-1:calcRegion(end)+1,i) = smooth(inputDensityArray(1,calcRegion(1)-1:calcRegion(end)+1,i),smoothAmount);
        inputRadiusArrayS(1,calcRegion(1)-1:calcRegion(end)+1,i) = inputRadiusArray(1,calcRegion(1)-1:calcRegion(end)+1,i);
    end
    
    if(zeroOn)
        %zero
        for i=1:length(inputDensityArray(1,1,:))
            shiftByThis(i) = mean(inputDensityArrayS(1,calcRegion(end)-8:calcRegion(end)-4,i));
            %shiftByThis(i) = min(inputDensityArrayS(1,calcRegion,i));
            %inputDensityArrayS(1,:,i) = inputDensityArrayS(1,:,i) - shiftByThis(i);
            shiftByThisMin(i) = min(inputDensityArrayS(1,calcRegion,i));
            inputDensityArrayS(1,:,i) = inputDensityArrayS(1,:,i) - shiftByThis(i);
        end
    end
else
    inputDensityArrayS = inputDensityArray;
    inputRadiusArrayS = inputRadiusArray;
    
    if(zeroOn)
        for i=1:length(inputDensityArray(1,1,:))
            shiftByThis(i) = mean(inputDensityArrayS(1,calcRegion(end)-5:calcRegion(end),i));
            inputDensityArrayS(1,:,i) = inputDensityArrayS(1,:,i) - shiftByThis(i);
        end        
    end
end

%convert to potential (via LDA)
radiusVector = inputRadiusArrayS(1,:,i); %x vector
radiusVectorMeters = radiusVector.*kpixelLength; %convert x vector to m
potential = 0.5 .* massL6 .* omegar^2 .* radiusVectorMeters.^2;


%Scale x&y to H.O. dimensionless
%Density SI = m^-2
radProfilesDensityConv = inputDensityArrayS./(massL6*omegaz / hbar);
%Potential SI = kg * s^-2 * m^2 // Energy units = kg * m^2 * s^-2
radProfilesPotentialConv = potential./((hbar*omegaz)/2);
% plot(radProfilesPotentialConv,radProfilesDensityConv);
figure(98); plot(potential,inputDensityArrayS,'.'); grid on; title('SI units density vs potential');
figure(99); plot(radProfilesPotentialConv,radProfilesDensityConv,'.'); grid on; title('H.O. units density vs potential');

%specify fitrange:
%TWEAK THIS:
%cutoffLeft = 40; cutoffRight = 12; (750G)
cutoffLeft = 40; cutoffRight = 12;
fitrange = calcRegion(end)-cutoffLeft:calcRegion(end)-cutoffRight;
figure(100); plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.'); grid on; title('H.O. units density vs potential fitrange');


%GG Fits:
%[a b] = fitGGComp(inputDensityArrayS,potential);



%Virial loop (bisection)
virialRsquares = []; tempGuess = []; virialTreturned = []; virialMu0returend = [];
virialFitFuncs = [];
for i=1:60;
tempGuess(i) = 10e-9 + 1*i*(1e-9); %Guess range in nK
betaeb = eb/(kB*tempGuess(i));

b2 = huicoeffs(betaeb,1,0);
b3 = huicoeffs(betaeb,2,0);

[f1, gof, foutput] = virial2Fit(radProfilesDensityConv(fitrange),radProfilesPotentialConv(fitrange),[b2 b3],omegaz);
%virial2Fit(radProfilesDensityAdjusted(fitrange),radProfilesPotential(fitrange),[b2 b3]);

%figure(i);
%plot(f1); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.'); hold off;

virialRsquares(i) = getfield(gof,'rsquare');
virialTreturned(i) = getfield(f1,'a');
virialMu0returend(i) = getfield(f1,'b');
virialFitFuncs{i} = f1;
end

if(0)
j = 1;
thresh = 0.003;
while bicond > thresh;
    tempGuessBi(j) = 10e-9; 
    tempHigh = 100e-9; %high guess initial
    tempLow = 10e-9; %low guess initial
    
    betaeb = eb/(kB*tempGuessBi(j));
    
    b2 = huicoeffs(betaeb,1,0);
    b3 = huicoeffs(betaeb,2,0);
    
    [f1, gof, foutput] = virial2Fit(radProfilesDensityConv(fitrange),radProfilesPotentialConv(fitrange),[b2 b3],omegaz);
    

    virialRsquares(j) = getfield(gof,'rsquare');
    virialTreturned(j) = getfield(f1,'a');
    virialMu0returend(j) = getfield(f1,'b');
    virialFitFuncs{j} = f1;
    
    %What of the returned temperature versus the guess?
    bicond = virialTreturned(j) - tempGuessBi(j);
    
    if(bicond < thresh)
        %congrats the returned is equal to the guess
        break;
    elseif(bicond < 0)
        %returned temperature larger than guess
        
    
    
    j = j+1;
    end
end

end



if(0)
    close all;
    showCurve = 30;
    plot(tempGuess,'.'); hold on; plot(virialTreturned,'.r');
    figure(2); plot(virialMu0returend.*(((hbar*omegaz)/2)),'.');
    figure(3); plot(virialRsquares,'.');    
    figure(4); plot(virialFitFuncs{showCurve}); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.');
    figure(5); plot(radProfilesPotentialConv,radProfilesDensityConv,'.'); hold on; plot(virialFitFuncs{showCurve}); grid on; axis([0 7 0 0.2]);
    plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'go');
end

TAbs = 17.5e-9;
muCloud = 1.798e-30;
%T / TF calc:
for i=1:length(inputDensityArray(1,1,:))
    densityCloud(1,:,i) = radProfilesDensityConv(1,:,i).*(massL6*omegaz / hbar)/2;
    for j=1:length(radProfilesDensityConv(1,:,1))
        EFermiAcrossCloud(j,i) = ((pi*hbar^2)/(massL6)).*densityCloud(1,j,i);
        TFermiAcrossCloud(j,i) = EFermiAcrossCloud(j,i)/kB;
        TonTFAcrossCloud(j,i) = TAbs / TFermiAcrossCloud(j,i);
        kfs(j,i) = real(sqrt(2*pi*densityCloud(1,j,i)));
        logkfa2ds(j,i) = log(kfs(j,i).*a2d);       
    end
end



lambdaDB = sqrt((2*pi*hbar^2)/(massL6*kB*TAbs));
%Density vs radius(meters):
%figure(222); plot(radiusVectorMeters,inputDensityArrayS(1,:,1),'.');


%Ideal chemical potential from atom number and given T: 
idealDensityArray = (2/(lambdaDB.^2)).*log(1 + exp((muCloud - potential)./(kB*TAbs)));
%plot(potential,idealDensityArray,'.'); hold on; grid on;
%plot(potential,inputDensityArrayS,'.r');

betaMuComplex = log(-1 + exp((lambdaDB^2) .* inputDensityArrayS(1,:,1)));
betaMu = real(betaMuComplex);
betaMuflip = fliplr(betaMu);

nOnn0 = (lambdaDB^2).*((inputDensityArrayS(1,:,1))./(log(1+exp(betaMu))));
nOnn0 = 2.*((log(1+exp(betaMu)))./(lambdaDB^2));
bulkn = log(1+exp(muCloud/(kB*TAbs)));
%figure(223); plot(betaMuflip,fliplr(nOnn0),'.');
%figure(225);
%plot(betaMuflip,fliplr(((lambdaDB^2).*((inputDensityArrayS(1,:,1))./bulkn))./2),'.');
if(0)
    figure(11); plot(betaMuflip,fliplr(nOnn0),'.');
    grid on;
    axis([-3 6 0 4]);
end


if(0)
    close all;
    %plot(TonTFAcrossCloud(1:30,1));
    fa = figure(133);
    legendString = [num2str(field) 'G   T_{abs} = ' num2str(TAbs)];
    plot(logkfa2ds(:,1),TonTFAcrossCloud(:,1),'DisplayName',legendString); grid on; axis([0 4 0 2]); title('T/TF vs log(kF a2d)');
    legend('show');
    
    figname = [num2str(field) '_TonTF_vs_logkFa2d_VirialFit'];
    figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\';
    
    saveas(fa,[figdirectory figname '.fig'],'fig');
    saveas(fa,[figdirectory figname '.png'],'png');

end

%EFermi = (hbar^2 * 4 * pi * n) / (2*massL6);
%TFermi = EFermi/kB;

%TonTF = 20.24e-9 / TFermi;






yarray = radProfilesDensityConv;
xarray = radProfilesPotentialConv;






end