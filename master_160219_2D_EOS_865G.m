directory = 'C:\Data\160219_2D_EOS_865G_approx10katoms_HighIntensity_freq_5p78kHz\';
date = '160219';
camera = 'topcam';
varstring = 'Isat';
%varstring2 = 'Holdtime';
pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
Isat = 134*10; %134*x us
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
widthYlogfile = []; widthXlogfile = []; nROIlogfile = [];
[fileLocList,varData] = generateFromLogfile(directory,date,varstring,camera);
widthYlogfile = varData(:,4);
widthXlogfile = varData(:,3);
nROIlogfile = varData(:,5);

imageArray = []; imageState1Array = []; beamImageArray = [];
%Pull images:
for i=1:length(fileLocList)
    Isat = varData(i,1);
    [beamImage,atom1Image,atom2Image] = PullSPE(fileLocList{i},Isat,OD);
    imageArray(:,:,i) = atom2Image(:,:,1);
    imageState1Array(:,:,i) = atom1Image(:,:,1);
    beamImageArray(:,:,i) = beamImage(:,:,1);    
end

%Crop images:
imageArrayC = []; imageArrayTC = [];
imageArrayHighInt = []; imageArrayLowInt = [];
ROIx = 1:length(imageArray(1,:,1));
ROIy = 1:length(imageArray(:,1,1));
CrossROIy = 30:160; 
CrossROIx = 40:180; %The cross is inside the region specified above.
CrossROIy = 1:100;
CrossROIx = 1:110;

TightROIx = 60:170;
TightROIy = 50:150;
%Split into high and low intensity arrays
imageArrayC = imageArray(TightROIy,TightROIx,:);
imageArrayTC = imageArray(TightROIy,TightROIx,:);

imageArrayHighInt = imageArrayC(:,:,1:22);
imageArrayLowInt = imageArrayC(:,:,23:end);

%Display every X image:
if(0)
for i=1:length(imageArray(1,1,:))
    if(mod(i,2) == 0)       
        figure(i);
        imagesc(imageArray(:,:,i));        
    end
end
for i=1:length(imageArrayC(1,1,:))
    if(mod(i,1) == 0)       
        figure(i);
        imagesc(imageArrayTC(:,:,i));        
    end
end
end

%Radially averaged profiles:
radProfiles = []; radProfilesT = []; center = [];
disp('Radially averaging...');
for i=1:length(imageArrayC(1,1,:))
    [radProfilesT(:,:,i),center(:,i)] = radAverageBigSquare(imageArrayC(:,:,i));
    radProfiles(:,:,i) = radProfilesT(:,1:end-5,i);
end

%Shift wings of radial profiles to zero:
shiftBy = [];
for i=1:length(radProfiles(1,1,:))
    shiftBy(i) = mean(radProfiles(1,72:84,i));
    radProfiles(1,:,i) = radProfiles(1,:,i) - shiftBy(i);
end

if(0)
    for i=1:length(imageArrayC(1,1,:))
        if(mod(i,2) == 0)
            figure(i);
            plot(fgr(gcoefsR(:,i),1:180),'g'); hold on;  plot(radProfiles(2,:,i),radProfiles(1,:,i),'r'); hold off;
        end
    end
end

%%%%%Fit Functions:
%gaussian with 4 variables:
fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4));
%gaussian with 2 variables (for radial profiles)
fgr = @(p,x)(p(1).*exp((-1).*((x).^2) ./ (2.*p(2).^2)));
%polylog order 1 function:
%fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*x.^2)./(p(3)^2))));
fgp = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*p(3).*x.^2)./(p(4)))));
%gaussian with 2 variables fixed 2.5 exponent:
fgr2p5 = @(p,x)(p(1).*exp((-1).*((x).^(2.5)) ./ (2.*p(2).^(2.5))));
%gaussian with 2 variables and fitted exponent:
fgrfe = @(p,x)(p(1).*exp((-1).*((x).^(p(2))) ./ (2.*p(3).^(2.5)))); 

gcoefsX = []; gcoefsY = []; centers = []; gcoefsXi = []; gcoefsYi = [];
gcoefsR = []; gcoefsR2p5 = []; sigmaR2p5 = []; gcoefsPolyLog1 = [];
gcoefsRFreeExp = []; sigmaR = [];
sigmaX = []; sigmaY = []; shiftFactor = []; shiftFactorR = [];
disp('Function Fitting...');
for i=1:length(imageArrayC(1,1,:))
    %Initial Fit for zeroing:
    gcoefsXi(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    gcoefsYi(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    gcoefsR(:,i) = gausFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    %gcoefsR2p5(:,i) = gausExp2p5FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    gcoefsPolyLog1(:,i) = polyLog1FitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i),'topcam');
    %gcoefsRFreeExp(:,i) = gausFreeExpFitHalf1D(radProfiles(1,:,i),radProfiles(2,:,i));
    
    %Shift wings to zero:
    shiftFactor(i) = (gcoefsXi(4,i)+gcoefsYi(4,i))/2;
    %imageArrayC(:,:,i) = imageArrayC(:,:,i) - shiftFactor(i);
        
    %Refit:
    gcoefsX(:,i) = gausFit1D(mean(imageArrayC(CrossROIy,:,i),1)); %mean averages over y
    %Profile: plot(mean(imageArrayC(:,:,i),1))
    gcoefsY(:,i) = gausFit1D(mean(imageArrayC(:,CrossROIx,i),2)); %mean averages over x
    %Profile: plot(mean(imageArrayC(:,:,i),2))
    
    centers(:,i) = [ceil(gcoefsX(2,i)), ceil(gcoefsY(2,i))]; %center = [x y]
    sigmaX(:,i) = gcoefsX(3,i);
    sigmaY(:,i) = gcoefsY(3,i);
    sigmaR(:,i) = gcoefsR(2,i);
    %sigmaR2p5(:,i) = gcoefsR2p5(2,i);
end


%show function fits:
if(0)
    
    %plot(fgp(gcoefsPolyLog1(:,i),1:180));
    
    %data fit:
    plot(varData(1:30,1),sigmaX(1:30),'.'); hold on; plot(varData(31:end,1),sigmaX(31:end),'.r');
end


%%%%%Atom numbers:
%Tight ROI array:
pixelCounts = [];
for i=1:length(imageArrayC(1,1,:))
    pixelCounts(i) = sum(sum(imageArrayC(:,:,i)));
end

%%%%%Temperatures:
TonTFs = []; Ts = [];
for i=1:length(imageArrayC(1,1,:))
    TonTFs(i) = 1/(log(1+exp(gcoefsPolyLog1(2,i)/gcoefsPolyLog1(3,i)^2)));
    Ts(i) = gcoefsPolyLog1(4,i);
end


%Re-order each array:
sortedVarData = []; indexs = []; sigmaRSort = [];
sigmaXSort = []; sigmaYSort = []; sigmaR2p5Sort = [];
sigmaRsmSort = []; TonTFsSort = []; imageArrayCSort = [];

reorder = 0;
if(reorder)
    [sortedVarData,indexs] = sort(varData);
%indexs(:,1) is a vector of the sort (varstring).
    for i=1:length(sigmaX)
        sigmaXSort(i) = sigmaX(indexs(i));
        sigmaYSort(i) = sigmaY(indexs(i));
        pixelCountsSort(i) = pixelCounts(indexs(i));
        sigmaRSort(i) = sigmaR(indexs(i));
        %sigmaR2p5Sort(i) = sigmaR2p5(indexs(i));
        %sigmaRsmSort(i) = sigmaRsm(indexs(i));
        TonTFsSort(i) = TonTFs(indexs(i));
        imageArrayCSort(:,:,i) = imageArrayC(:,:,indexs(i));
    end
else
    sortedVarData = varData;
    for i=1:length(sigmaX)
        sigmaXSort(i) = sigmaX(i);
        sigmaYSort(i) = sigmaY(i);
        pixelCountsSort(i) = pixelCounts(i);
        sigmaRSort(i) = sigmaR(i);
        %sigmaR2p5Sort(i) = sigmaR2p5(indexs(i));
        %sigmaRsmSort(i) = sigmaRsm(indexs(i));
        TonTFsSort(i) = TonTFs(i);
        imageArrayCSort(:,:,i) = imageArrayC(:,:,i);
    end
end


%Average over same motfet data points:
%Bug fix history: 150109 fixed endNum
j=1; runTotal = 0; motFets = []; widthsX = []; widthsY = [];
stdDevWidthsX = []; stdDevWidthsY = []; pixelNumbers = [];
pixelNumbersStdDev = []; sMomentX = []; sMomentY = [];
widthsR = []; stdDevWidthsR = []; widthsR2p5 = [];
widthsPsm = []; stdDevWidthsPsm = []; TonTFsm = []; stdDevTonTFsm = [];
sMomentXStdDev = []; sMomentYStdDev = []; stdDevWidthsR2p5 = [];
imageArrayAvgsNoCenter = [];
prev = sortedVarData(1); imageArrayAvgs = []; runTotals = [];
for i=1:length(sortedVarData)
    curr = sortedVarData(i);
    
    if( curr == prev )
        runTotal = runTotal+1;
    else
        %hit next value
        runTotals(j) = runTotal;
        startNum = i-runTotal;
        endNum = i-1;
        widthsX(j) = mean(sigmaXSort(i-runTotal:endNum));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:endNum));
        widthsY(j) = mean(sigmaYSort(i-runTotal:endNum));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:endNum));
        widthsR(j) = mean(sigmaRSort(i-runTotal:endNum));
        stdDevWidthsR(j) = std(sigmaRSort(i-runTotal:endNum));
        %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
        %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
        %widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        %stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:endNum));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:endNum));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:endNum));
        imageArrayAvgsNoCenter(:,:,j) = mean(imageArrayCSort(:,:,i-runTotal:endNum),3);
        
        motFets(j) = sortedVarData(i-runTotal);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:endNum));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:endNum));
        runTotal = 1;
        j = j+1;
    end
    if( i == length(sortedVarData))
        %Final run:
        if(i == runTotal)
            runTotal = runTotal-1;
        end
        runTotals(j) = runTotal;
        widthsX(j) = mean(sigmaXSort(i-runTotal:i));
        stdDevWidthsX(j) = std(sigmaXSort(i-runTotal:i));
        widthsY(j) = mean(sigmaYSort(i-runTotal:i));
        stdDevWidthsY(j) = std(sigmaYSort(i-runTotal:i));
        widthsR(j) = mean(sigmaRSort(i-runTotal:i));
        stdDevWidthsR(j) = std(sigmaRSort(i-runTotal:i));
        %widthsR2p5(j) = mean(sigmaR2p5Sort(i-runTotal:i));
        %stdDevWidthsR2p5(j) = std(sigmaR2p5Sort(i-runTotal:i));
        %widthsPsm(j) = mean(sigmaRsmSort(i-runTotal:i));
        %stdDevWidthsPsm(j) = std(sigmaRsmSort(i-runTotal:i));
        
        TonTFsm(j) = mean(TonTFsSort(i-runTotal:i));
        stdDevTonTFsm(j) = std(TonTFsSort(i-runTotal:i));
        
        imageArrayAvgs(:,:,j) = centerAndAverage(imageArrayCSort(:,:,i-runTotal:i));
        imageArrayAvgsNoCenter(:,:,j) = mean(imageArrayCSort(:,:,i-runTotal:i),3);
        
        motFets(j) = sortedVarData(i);
        pixelNumbers(j) = mean(pixelCountsSort(i-runTotal:i));
        pixelNumbersStdDev(j) = std(pixelCountsSort(i-runTotal:i));
        runTotal = 1;
        j = j+1;
    end
    
    prev = curr;
end

%Radially averaged average profiles:
radProfilesAvg = []; radProfilesTAvg = []; centerAvgs = [];
disp('Radially averaging...');
for i=1:length(imageArrayAvgs(1,1,:))
    [radProfilesTAvg(:,:,i),centerAvgs(:,i)] = radAverageBigSquare(imageArrayAvgs(:,:,i));
    radProfilesAvg(:,:,i) = radProfilesTAvg(:,1:end-5,i);
end

%Shift wings of radial profiles to zero:
shiftByAvgs = [];
for i=1:length(radProfilesAvg(1,1,:))
    shiftByAvgs(i) = mean(radProfilesAvg(1,72:87,i));
    radProfilesAvg(1,:,i) = radProfilesAvg(1,:,i) - shiftByAvgs(i);
end

%Shift images by found factor:
imageArrayAvgsZerod = [];
for i=1:length(imageArrayAvgs(1,1,:))
    imageArrayAvgsZerod(:,:,i) = imageArrayAvgs(:,:,i) - shiftByAvgs(i);
end


if(0)
    for i=1:length(imageArrayAvgs(1,1,:))
        if(mod(i,1) == 0)
            figure(i);
            hold on; plot(radProfilesAvg(2,:,i),radProfilesAvg(1,:,i),'r'); line([0 100],[0 0]); hold off;
        end
    end
end

radproftoFityReal = radProfilesAvg(1,:,2)./(pixelLength^2);
%radproftoFityRealConv = radproftoFityReal.*10^(-11);
%radproftoFitxReal = radProfilesAvg(2,58:76,2).*pixelLength;

%Fit Virial:
massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
magfield = 865;
omegaz = 5789 * 2 * pi;
omegar = 25 * 2 * pi; %Double check this weak trapping
az = sqrt(hbar/(massL6*omegaz));
pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
kpixelLength = (13e-6*(83/400));
radiusVector = 1:(length(radProfilesAvg(2,:,1)));
radiusVectorMeters = radiusVector.*kpixelLength;
potential = 0.5 .* massL6 .* omegar^2 .* radiusVectorMeters.^2;
potential_nk = (potential./kB).*(10^9);

%convert radProfile pixels to potential:
radProfilesPotential = 0.5 .* massL6 .* omegar^2 .*(radProfilesAvg(2,:,2).*pixelLength).^2;
%plot(radProfilesPotential,radProfilesAvg(1,:,2)./(pixelLength^2));
%radProfilesPotentialConv = radProfilesPotential.*10^30;
radProfilesDensityAdjusted = radproftoFityReal.*2.*1.2; %2 for spin states and 1.2 for correction factor

a0 = 5.29e-11; %o.0
a2d = 47186.7 * a0; %From mathematica spreadsheet TODO generate text file
eb = (hbar^2)/(a2d^2 * massL6);


%Scale x&y to H.O. dimensionless 
%Density SI = m^-2
radProfilesDensityConv = radProfilesDensityAdjusted./(massL6*omegaz / hbar);
%Potential SI = kg * s^-2 * m^2 // Energy units = kg * m^2 * s^-2
radProfilesPotentialConv = radProfilesPotential./((hbar*omegaz)/2); 

fitrange = 60:75;

virialRsquares = []; tempGuess = []; virialTreturned = []; virialMu0returend = [];
virialFitFuncs = [];
for i=1:60;
tempGuess(i) = 10e-9 + i*(0.5e-9); %Guess 20nK
betaeb = eb/(kB*tempGuess(i));

b2 = huicoeffs(betaeb,1,0);
b3 = huicoeffs(betaeb,2,0);

[f1, gof, foutput] = virial2Fit(radProfilesDensityConv(fitrange),radProfilesPotentialConv(fitrange),[b2 b3]);
%virial2Fit(radProfilesDensityAdjusted(fitrange),radProfilesPotential(fitrange),[b2 b3]);

%figure(i);
%plot(f1); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.'); hold off;

virialRsquares(i) = getfield(gof,'rsquare');
virialTreturned(i) = getfield(f1,'a');
virialMu0returend(i) = getfield(f1,'b');
virialFitFuncs{i} = f1;
end

if(0)
    close all;
    plot(tempGuess,'.'); hold on; plot(virialTreturned,'.r');
    figure(2); plot(virialRsquares,'.');
    
    figure(3); plot(virialFitFuncs{19}); hold on; plot(radProfilesPotentialConv(fitrange),radProfilesDensityConv(fitrange),'.');
end

%T / TF calc:

n = 0.25*(massL6*omegaz / hbar)/2;
EFermi = (hbar^2 * 4 * pi * n) / (2*massL6);
TFermi = EFermi/kB;

TonTF = 20.24e-9 / TFermi;

atomNum = [];
for i=1:length(imageArrayAvgs(1,1,:))
    atomNum(i) = sum(sum(imageArrayAvgs(:,:,i)));
end

%{
%----------------------------------------------------------------%
%virial_fit_residuals(0.2, 0.6, eb, density(cutoff:endp), potential(cutoff:endp));


    %Define BetaEb guess range:
    beta_eb = [];
    beta_eb(1) = beta_eb_lower;
    beta_eb(2) = beta_eb_upper;
    beta_eb(3) = beta_eb(1) + ((beta_eb_upper - beta_eb_lower) / 2);
   
    b2 = []; %Array of size 3
    b3 = []; %Array of size 3

       % Use virial coeffcients to find T and mu0 with virial coeffcients
    def v_residuals(params, x, b2, b3, eps_data, y=None):
        %T = params['T'].value;
        %mu0 = params['mu0'].value;
        %alpha = params['alpha'].value
        %dens_fit = ((2/((2*pi*hbar**2)/((mass_li*kb*T)))) * 
        %           (np.exp((1/(kb*T)) * (mu0 - x))) + 
        %           (2*b2*np.exp(2*(1/(kb*T))*(mu0 -x))) + 
        %           (3*b3*np.exp(3*(1/(kb*T))*(mu0 - x))))

        dens_fit = ((2/((2*pi*hbar^2)/(mass_li*kb*T))) * ...
                   (np.log( 1 + np.exp((1/(kb*T)) * (mu0 - x))) + ...
                   (2*b2*np.exp(2*(1/(kb*T))*(mu0 - x))) + ...
                   (3*b3*np.exp(3*(1/(kb*T))*(mu0 - x)))));
        if y is None:
            %return dens_fit
        %return (y - dens_fit) 

    eps_data = 1e-4;
    % Repeat bisection procedure for 10 iterations
    for i=1:25
        %i is the enumeration, beb the value; i=1:3
        for j=1:length(beta_eb) 
           %run virial code to get b2 and b3 coeff
           b2_proc = huicoeff(beta_eb(j),1,0); 
           b3_proc = huicoeff(beta_eb(j),2,0);

           %b2_out = b2_proc.communicate()[0]
           %b3_out = b3_proc.communicate()[0]

           b2(j) = b2_proc;
           b3(j) = b3_proc;
        end
        beta_eb_fit = unp.uarray(np.zeros(3),np.zeros(3))
        beta_eb_diff = unp.uarray(np.zeros(3),np.zeros(3))
        % fit for upper, mid and lower beta_eb
        for k=1:length(beta_eb)
            params = lmfit.Parameters() %lmfit is fitting routine
            params.add('T', value = 50e-9, min=5e-9, max=70e-9)
            params.add('mu0', value = 1.7e-30, min=1.e-30, max=4e-30)
            %params.add('alpha', value = 1.25, min = 1, max=2)
            fit_output = lmfit.minimize(v_residuals, params, 
                                 args=(potential, b2[k], b3[k], eps_data, density))

            beta_fit = 1/(kb * params['T'].value)
            beta_eb_fit[k] = eb * beta_fit
            beta_eb_diff[k] = np.abs(beta_eb_fit[k] - beb)
            %print(lmfit.fit_report(fit_output))

            %print(beta_eb_diff[i]*100)

        if (beta_eb_diff[0] - beta_eb_diff[1]) < (beta_eb_diff[2] -
                beta_eb_diff[1]):
            beta_eb[2] = beta_eb[1]
            beta_eb[1] = beta_eb[0] + ((beta_eb[2] - beta_eb[0]) / 2)
        else:
            beta_eb[0] = beta_eb[1]
            beta_eb[1] = beta_eb[0] + ((beta_eb[2] - beta_eb[0]) / 2)
        end
        end
    end
        
    beta_eb_max = eb * 1/(kb * (params['T'].value+params['T'].stderr))
    beta_eb_err = np.abs(beta_eb_max - beta_eb[1])
    beta_eb_final = eb * 1/(kb * (params['T'].value))
    fit_output = lmfit.minimize(v_residuals, params, 
                                args=(potential, b2[1], b3[1], eps_data, density))
    
    print('Final betaEb: {0:f} +/- {1:f}'.format(beta_eb_final, beta_eb_err))
    print(lmfit.fit_report(fit_output))
    fit = v_residuals(params, potential, b2[1], b3[1], 1)
    return fit, fit_output%, beta_eb[1]  
        end

%}