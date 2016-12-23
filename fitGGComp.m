%function [xarray, yarray] = fitGGComp(inputDensityArray, inputPotentialArray)
%Take input density vs potential experimental data and fit to a range of
%functions of f_n vs beta_mu which are virial+GG composite functions

%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\fitFuncsCompCell.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\fitFuncsCompCell_betaMuMax1.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\fitFuncsCompCell_161129.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\fitFuncsCompCell_151216.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\betaEb_virialsArray.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\betaEb_virialsArray_151216.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\betamuggs.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\fnggs.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\fnggs_151216.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\betamuggs_151216.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\kdata_880_density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\kdata_880_potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\kdata_865_density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\kdata_865_potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\kdata_920_density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\kdata_920_potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_160920_690G_7k_Density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_160920_690G_7k_Potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161003_725G_10k_Density.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161003_725G_10k_Potential.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161007_750G_10k_Density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161007_750G_10k_Potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161006_712G_10k_Density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161006_712G_10k_Potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161006_706G_10k_Density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161006_706G_10k_Potential_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161004_735G_10k_Density_SI.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\Density_Profiles\ColOsc_161004_735G_10k_Potential_SI.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\virialtwos.mat');
%load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\virialthrees.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\virialtwos_151216.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\virialthrees_151216.mat');
load('C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\composite_curve_data\GGSurface.mat');



massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1

close all;
%choices = [1,2,3,5,6,7,9]
for choices = [8]
%--------------Init variables:
%<<<<<<<
datasetChoice = choices; % 1 = 880G, 2 = 920G, 3 = 865G
%betaEbEval = 5; %13 = 0.26, 23 = 0.46, 3 = 0.06
%zeroCorrect = 1; %re-zero data?
%<<<<<<<

%--------------Dataset creation (from *.mat loaded density/potentials):
dataset2D{1}.name = '880 Kristian EOS'; %Kristian result: 17.5nK +-4
dataset2D{1}.a2d = 3.9910e-06;
dataset2D{1}.density = kdata_880_density_SI;
dataset2D{1}.potential = kdata_880_potential_SI;
dataset2D{1}.fitpoints = 38:60;
dataset2D{1}.rezero.yesno = 0;
dataset2D{1}.rezero.endminusleft = 1;
dataset2D{1}.rezero.endminusright = 1;
dataset2D{1}.betaEBidx = 14;

dataset2D{2}.name = '920 Kristian EOS'; %Kristian result: 18.5nK +-3
dataset2D{2}.a2d = 8.9761e-06;
dataset2D{2}.density = kdata_920_density_SI;
dataset2D{2}.potential = kdata_920_potential_SI;
dataset2D{2}.fitpoints = 36:50;
dataset2D{2}.rezero.yesno = 0;
dataset2D{2}.rezero.endminusleft = 0;
dataset2D{2}.rezero.endminusright = 0;
dataset2D{2}.betaEBidx = 3;

dataset2D{3}.name = '865 Kristian EOS'; %Kristian result: 22nK +-4
dataset2D{3}.a2d = 2.7867e-06;
dataset2D{3}.density = kdata_865_density_SI;
dataset2D{3}.potential = kdata_865_potential_SI;
dataset2D{3}.fitpoints = 32:59;
dataset2D{3}.rezero.yesno = 0;
dataset2D{3}.rezero.endminusleft = 0;
dataset2D{3}.rezero.endminusright = 0;
dataset2D{3}.betaEBidx = 23;

dataset2D{4}.name = 'ColOsc 160920 690G 7k'; %Resonance OD corr of Coll Osc
dataset2D{4}.a2d = 1.1220e-06;
dataset2D{4}.density = ColOsc_160920_690G_7k_Density;
dataset2D{4}.potential = ColOsc_160920_690G_7k_Potential;
dataset2D{4}.fitpoints = 27:33;
dataset2D{4}.rezero.yesno = 0;
dataset2D{4}.rezero.endminusleft = 0;
dataset2D{4}.rezero.endminusright = 0;
dataset2D{4}.betaEBidx = 38; %NO BETAEB DATA FOR THIS SET 

dataset2D{5}.name = 'ColOsc 161003 725G 10k'; %725G OD corr of Coll Osc
dataset2D{5}.a2d = 5.6412e-06;
dataset2D{5}.density = ColOsc_161003_725G_10k_Density;
dataset2D{5}.potential = ColOsc_161003_725G_10k_Potential;
dataset2D{5}.fitpoints = 32:46;
dataset2D{5}.rezero.yesno = 1;
dataset2D{5}.rezero.endminusleft = 40;
dataset2D{5}.rezero.endminusright = 30;
dataset2D{5}.betaEBidx = 5;

dataset2D{6}.name = 'ColOsc 161007 750G 10k'; %750G OD corr of Coll Osc
dataset2D{6}.a2d = 1.2421e-05;
dataset2D{6}.density = ColOsc_161007_750G_10k_Density_SI;
dataset2D{6}.potential = ColOsc_161007_750G_10k_Potential_SI;
dataset2D{6}.fitpoints = 35:46;
dataset2D{6}.rezero.yesno = 1;
dataset2D{6}.rezero.endminusleft = 30;
dataset2D{6}.rezero.endminusright = 21;
dataset2D{6}.betaEBidx = 1;

dataset2D{7}.name = 'ColOsc 161006 712G 10k'; %712G OD corr of Coll Osc
dataset2D{7}.a2d = 3.3535e-06;
dataset2D{7}.density = ColOsc_161006_712G_10k_Density_SI;
dataset2D{7}.potential = ColOsc_161006_712G_10k_Potential_SI;
dataset2D{7}.fitpoints = 34:48;
dataset2D{7}.rezero.yesno = 1;
dataset2D{7}.rezero.endminusleft = 30;
dataset2D{7}.rezero.endminusright = 21;
dataset2D{7}.betaEBidx = 13;

dataset2D{8}.name = 'ColOsc 161006 706G 10k'; %706G OD corr of Coll Osc
dataset2D{8}.a2d = 2.5560e-06;
dataset2D{8}.density = ColOsc_161006_706G_10k_Density_SI;
dataset2D{8}.potential = ColOsc_161006_706G_10k_Potential_SI;
dataset2D{8}.fitpoints = 32:42;
dataset2D{8}.rezero.yesno = 1;
dataset2D{8}.rezero.endminusleft = 37;
dataset2D{8}.rezero.endminusright = 30;
dataset2D{8}.betaEBidx = 28; %NO BETAEB DATA FOR THIS SET 

dataset2D{9}.name = 'ColOsc 161004 735G 10k'; %735G OD corr of Coll Osc
dataset2D{9}.a2d = 7.9622e-06;
dataset2D{9}.density = ColOsc_161004_735G_10k_Density_SI;
dataset2D{9}.potential = ColOsc_161004_735G_10k_Potential_SI;
dataset2D{9}.fitpoints = 30:39;
dataset2D{9}.rezero.yesno = 0;
dataset2D{9}.rezero.endminusleft = 37;
dataset2D{9}.rezero.endminusright = 30;
dataset2D{9}.betaEBidx = 2;

%--------------Surface Init
surf(GGSurface{1}.x,GGSurface{1}.y,GGSurface{1}.z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([-5 2 0 1 0 1.5]);

surface_betaebs = GGSurface{1}.y(:,1);
surface_betamus = GGSurface{1}.x(1,:);
surface_fns = GGSurface{1}.z;


%--------------Pre-computation:
betaEbEval = dataset2D{datasetChoice}.betaEBidx;

fitIndexsBetaMu = [];
fitIndexsBetaMu = dataset2D{datasetChoice}.fitpoints;
eb = (hbar^2)/((dataset2D{datasetChoice}.a2d^2) * massL6);
%(betaebvirials(betaEbEval,1))
thisBetaEbTemperature = (1/(kB*((betaebvirials(betaEbEval,1))/eb)));

compxs = -10:0.1:1;
virialxs = -10:0.01:5;

%inputDensityArray = kdata_880_density_SI;
%inputPotentialArray = kdata_880_potential_SI;
inputDensityArray = dataset2D{datasetChoice}.density;
inputPotentialArray = dataset2D{datasetChoice}.potential;

figure(1);
plot(inputPotentialArray,inputDensityArray,'.');
grid on;

if(dataset2D{datasetChoice}.rezero.yesno)
    shiftByThis = []; shiftByThisMin = [];
    for i=1:length(inputDensityArray)
        shiftByThis(i) = mean(inputDensityArray(end-dataset2D{datasetChoice}.rezero.endminusleft:end-dataset2D{datasetChoice}.rezero.endminusright));
        %shiftByThisMin(i) = min(inputDensityArray);
        inputDensityArray(i) = inputDensityArray(i) - shiftByThis(i);
    end   
end
 
%xdataorig = (guessMu_O - inputPotentialArray).*(beta);
%ydataorig = (inputDensityArray.*(lambdaDB.^2))./(2*pi);

%----------------Gain equation GG+virial composite in a string:

eq = formula(fitFuncsComp{betaEbEval}); %Formula of fitted equation
parameters = coeffnames(fitFuncsComp{betaEbEval}); %All the parameter names
values = coeffvalues(fitFuncsComp{betaEbEval}); %All the parameter values
for idx = 1:numel(parameters)
      param = parameters{idx};
      l = length(param);
      loc = regexp(eq, param); %Location of the parameter within the string
      while ~isempty(loc)     
          %Substitute parameter value
          eq = [eq(1:loc-1) num2str(values(idx)) eq(loc+l:end)];
          loc = regexp(eq, param);
          if param == 'e'
              loc = []; %the second 'e' is in 'exp' /facepalm
          end
      end
end

%ydata = feval(fitFuncsComp{betaEbEval},inputPotentialArray);


%----------------Main fitting routine:
%x = 1,2,3,...,length(inputDensityArray)
%feval(fitFuncsComp{betaEbEval},(mu0guess - inputPotentialArray(X)).*betaguess)

minimizefun = @(p,x)((inputDensityArray(x).*((hbar^2)./(massL6)).*(1/(kB.*p(1).*1e-9)))' - feval(fitFuncsComp{betaEbEval},(p(2).*(kB*1e-9) - inputPotentialArray(x)).*(1/(kB.*p(1).*1e-9))));

mu0guessConv = 150;
betaguessConv = 35;

xConvFull = []; yConvFull = [];
xConvFull = (mu0guessConv.*(kB*1e-9) - inputPotentialArray).*(1/(kB.*betaguessConv.*1e-9));
yConvFull = (inputDensityArray.*(sqrt(((1/(kB.*betaguessConv.*1e-9)).*2*pi*hbar^2)/(massL6)).^2)./(2*pi));
%Find beta_mu's over which to fit:
%fitIndexsBetaMu = find(yConvFull < feval(fitFuncsComp{betaEbEval},1) & yConvFull > feval(fitFuncsComp{betaEbEval},-4));
%fitIndexsBetaMu = find(xConvFull < 1.5 & xConvFull > -3);


xdatafit = fitIndexsBetaMu; %fit on these beta mu points
ydatafit = zeros(length(xdatafit),1); %we want to fit to zeros
p0 = [betaguessConv mu0guessConv]; %[betaguess mu0guess]
lb = [1 1];
ub = [30 1000];
curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-20,'TolX',1e-20);
%curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',100000);
[f1,gof,r,exitflag,output,lambda,J] = lsqcurvefit(minimizefun,p0,xdatafit',ydatafit,lb,ub,curvefitoptions);
outputCoefError = nlparci(f1,r,'jacobian',J);

xwithfit = (f1(2).*(kB*1e-9) - inputPotentialArray).*(1/(kB.*f1(1).*1e-9));
ywithfit = (inputDensityArray.*(sqrt(((1/(kB.*f1(1).*1e-9))*2*pi*hbar^2)/(massL6)).^2)./(2*pi));

fittedTemperature = f1(1);
fittedMu0 = f1(2);


%----------------Calculation of values requiring T & mu0
TAbs = fittedTemperature * 10^(-9);
muCloud = fittedMu0;
%T / TF calc:
EFermiAcrossCloud = []; TFermiAcrossCloud = []; TonTFAcrossCloud = [];
kfs = []; logkfa2ds = []; densityCloud = [];
for i=1
    densityCloud(1,:,i) = inputDensityArray;
    for j=1:length(inputDensityArray)
        EFermiAcrossCloud(j,i) = ((pi*hbar^2)/(massL6)).*densityCloud(1,j,i);
        TFermiAcrossCloud(j,i) = EFermiAcrossCloud(j,i)/kB;
        TonTFAcrossCloud(j,i) = TAbs / TFermiAcrossCloud(j,i);
        kfs(j,i) = real(sqrt(2*pi*densityCloud(1,j,i)));
        logkfa2ds(j,i) = log(kfs(j,i).*a2d);       
    end
end

%-----------------Density/Density0 vs BetaMu:

Fn = []; BetaMu = []; denomNN0 = []; NonN0 = []; 

Fn = ywithfit(fitIndexsBetaMu) %Need Fp' not Fn...
BetaMu = xwithfit(fitIndexsBetaMu);

virialx = virialxs;
virialy = virialthrees(betaEbEval,:);

ggy = feval(fitFuncsComp{betaEbEval},xwithfit);
ggx = xwithfit;

%Fn = ggy; %Need Fp' not Fn...
%BetaMu = ggx;

denomNN0 = 2.*log(1+exp(BetaMu));

for i=1:length(Fn)
    denomNN0(i) = 2.*log(1+exp(BetaMu(i)));
    NonN0(i) = (Fn(i)./denomNN0(i)).*4.*pi;
end

fignumberNN0 = 1212 + choices;
figure(fignumberNN0);
plot(BetaMu,NonN0,'.');
grid on;
axis([-6 4 0 4]);

%figure(555);
%plot(inputPotentialArray(1:30),TonTFAcrossCloud(1:30),'.');
%title('T on T_F vs potential across the cloud');

%---------------Final Graph:


%close all;
fignumber = 222 + choices;
figure(fignumber);
plot(xwithfit,ywithfit,'.-');
hold on;
plot(xwithfit,feval(fitFuncsComp{betaEbEval},xwithfit),'-');
plot(xwithfit(fitIndexsBetaMu),ywithfit(fitIndexsBetaMu),'og');
grid on;
plot(virialxs,virialtwos(betaEbEval,:),'--');
plot(virialxs,virialthrees(betaEbEval,:),'-.');
plot(betaMusGG{betaEbEval},fnGG{betaEbEval},'-.'); %real GG points
%titlestring{choices} = {[' ' dataset2D{datasetChoice}.name ' Fit       \beta\E_b = ' num2str(betaebvirials(betaEbEval,1))  ...
%        ' (T_{beb} = ' num2str(thisBetaEbTemperature./(1e-9)) ' nK)'],['T_{fit} = ' num2str(f1(1)) ' (' num2str(outputCoefError(1,1)) ' ' num2str(outputCoefError(1,2)) ') nK     ' ...
%        ' \mu_0 = ' num2str(f1(2)) ' (' num2str(outputCoefError(2,1)) ' ' num2str(outputCoefError(2,2)) ') nK'],['Fit quality = ' num2str(1-gof) '      T/T_F = ' num2str(mean(TonTFAcrossCloud(1:4,1))) ' ']};
titlestring{choices} = {[' ' dataset2D{datasetChoice}.name ' Fit']};
textstring{choices} = {'Properties: ',['\beta E_B = ' num2str(betaebvirials(betaEbEval,1),4) ' (T_{\beta E_B} = ' num2str(thisBetaEbTemperature./(1e-9),4) ' nK)'],['T_{fit} = ' num2str(f1(1),4) ' (' num2str(outputCoefError(1,1),4) ', ' num2str(outputCoefError(1,2),4) ') nK     ' ...
        ''],['\mu_0 = ' num2str(f1(2),4) ' (' num2str(outputCoefError(2,1),4) ', ' num2str(outputCoefError(2,2),4) ') nK'],['Fit quality = ' num2str(1-gof,4) ''],['T/T_F = ' num2str(mean(TonTFAcrossCloud(1:4,1)),4) ' ']};
title(titlestring{choices});
text(-5.5,0.5,textstring{choices});
xlabel('\beta\mu');
ylabel('F_N');
axis([-6 4 -0.08 1]);

end


if(0)
%----------------Brute force method (old):
signature4to3 = []; signature4to2 = []; signature4to1 = []; signature4to0 = [];
betaEbEval = 13;
figure(1111);
plot(compxs,feval(fitFuncsComp{10}, compxs),'.');
hold on;
plot(compxs,feval(fitFuncsComp{12}, compxs),'.');
plot(compxs,feval(fitFuncsComp{13}, compxs),'.');
k = 1;
for i=1:50
    for j=1:30
    
    mu0ThisGuess = 1.5e-30 + i*1e-32;
    %mu0ThisGuess = 1.5e-30;
    betaThisGuess = 3e30 + j*0.5e29;
    
    xThis = (mu0ThisGuess - inputPotentialArray).*(betaThisGuess);
    yThis = (inputDensityArray.*((hbar^2)/(massL6)).*betaThisGuess);
    
    plot(xThis,yThis);
    
    %calculate residuals (subtract data from fit) at each point -4 to -1
    %and -4 to -2 and -4 to 0 (betaMu)...
    %'fit signature'
    
    pointsToSig4to3 = find(xThis < -3 & xThis > -4);
    pointsToSig4to2 = find(xThis < -2 & xThis > -4);
    pointsToSig4to1 = find(xThis < -1 & xThis > -4);
    pointsToSig4to0 = find(xThis < -0 & xThis > -4);
    
    signature4to3{k}.sig = yThis(pointsToSig4to3) - feval(fitFuncsComp{betaEbEval},xThis(pointsToSig4to3))';
    signature4to3{k}.mu0 = mu0ThisGuess;
    signature4to3{k}.beta = betaThisGuess;
    signature4to3{k}.x = xThis(pointsToSig4to3);
    signature4to3{k}.y = yThis(pointsToSig4to3);
    signature4to3{k}.score = abs(1/sum(signature4to3{k}.sig(:)));
    surface4to3x(i) = mu0ThisGuess;
    surface4to3y(j) = betaThisGuess;
    surface4to3z(i,j) = signature4to3{k}.score;
    
    signature4to2{k}.sig = yThis(pointsToSig4to2) - feval(fitFuncsComp{betaEbEval},xThis(pointsToSig4to2))';
    signature4to2{k}.mu0 = mu0ThisGuess;
    signature4to2{k}.beta = betaThisGuess;
    signature4to2{k}.x = xThis(pointsToSig4to2);
    signature4to2{k}.y = yThis(pointsToSig4to2);
    signature4to2{k}.score = abs(1/sum(signature4to2{k}.sig(:))); 
    surface4to2x(i) = mu0ThisGuess;
    surface4to2y(j) = betaThisGuess;
    surface4to2z(i,j) = signature4to2{k}.score;
    
    signature4to1{k}.sig = yThis(pointsToSig4to1) - feval(fitFuncsComp{betaEbEval},xThis(pointsToSig4to1))';
    signature4to1{k}.mu0 = mu0ThisGuess;
    signature4to1{k}.beta = betaThisGuess;
    signature4to1{k}.x = xThis(pointsToSig4to1);
    signature4to1{k}.y = yThis(pointsToSig4to1);
    signature4to1{k}.score = abs(1/sum(signature4to1{k}.sig(:)));
    surface4to1x(i) = mu0ThisGuess;
    surface4to1y(j) = betaThisGuess;
    surface4to1z(i,j) = signature4to1{k}.score;
    
    signature4to0{k}.sig = yThis(pointsToSig4to0) - feval(fitFuncsComp{betaEbEval},xThis(pointsToSig4to0))';
    signature4to0{k}.mu0 = mu0ThisGuess;
    signature4to0{k}.beta = betaThisGuess;
    signature4to0{k}.x = xThis(pointsToSig4to0);
    signature4to0{k}.y = yThis(pointsToSig4to0);
    signature4to0{k}.score = abs(1/sum(signature4to0{k}.sig(:)));
    surface4to0x(i) = mu0ThisGuess;
    surface4to0y(j) = betaThisGuess;
    surface4to0z(i,j) = signature4to0{k}.score;    
    
    k = k+1;
    end
end
grid on;

%plot score
for i=1:length(signature4to2)
    scores4to3(i) = signature4to3{i}.score;
    scores4to2(i) = signature4to2{i}.score;
    scores4to1(i) = signature4to1{i}.score;
    scores4to0(i) = signature4to0{i}.score;
end
figure(402);
plot(scores4to3,'.');
hold on;
plot(scores4to2,'.');
plot(scores4to1,'.');
plot(scores4to0,'.');
grid on;

%where is the score maximum?
maxScoreIndex = find(scores4to3 == max(scores4to3));
thisIndex = find(scores4to1 > 500);
goodIndex = 1219;

figure(13121);
plot(signature4to3{maxScoreIndex}.x,signature4to3{maxScoreIndex}.sig);
figure(13122);
%plot(signature4to3{maxScoreIndex}.x,signature4to3{maxScoreIndex}.y);
hold on;
plot(signature4to3{maxScoreIndex}.x,feval(fitFuncsComp{betaEbEval},signature4to3{maxScoreIndex}.x),'o');
plot((signature4to3{maxScoreIndex}.mu0 - inputPotentialArray).*(signature4to3{maxScoreIndex}.beta),(inputDensityArray.*((hbar^2)/(massL6)).*signature4to3{maxScoreIndex}.beta));
plot((signature4to3{maxScoreIndex}.mu0 - inputPotentialArray).*(signature4to3{maxScoreIndex}.beta),feval(fitFuncsComp{betaEbEval},(signature4to3{maxScoreIndex}.mu0 - inputPotentialArray).*(signature4to3{maxScoreIndex}.beta)));
grid on;


figure(13124);
%plot(signature4to3{maxScoreIndex}.x,signature4to3{maxScoreIndex}.y);
hold on;
plot(signature4to1{goodIndex}.x,signature4to1{goodIndex}.y,'o');
plot((signature4to1{goodIndex}.mu0 - inputPotentialArray).*(signature4to1{goodIndex}.beta),(inputDensityArray.*((hbar^2)/(massL6)).*signature4to1{goodIndex}.beta),'.');
plot((signature4to1{goodIndex}.mu0 - inputPotentialArray).*(signature4to1{goodIndex}.beta),feval(fitFuncsComp{betaEbEval},(signature4to1{goodIndex}.mu0 - inputPotentialArray).*(signature4to1{goodIndex}.beta)));
grid on;
title(['880G Fit \beta\E_b = ' num2str(betaebvirials(betaEbEval,1)) ',  T = ' num2str((1/signature4to1{goodIndex}.beta)/kB) 'K,  \mu_0 = ' num2str(signature4to1{goodIndex}.mu0) ' J, Virial from -3.5 \beta\mu']);
axis([-4.2 2 -0.08 0.5]);

%score surface, based on mu0 and beta...
figure(333);

surf(surface4to3y,surface4to3x,surface4to3z);
title('Fit from -4 to -3 betamu');
figure(334);
surf(surface4to2y,surface4to2x,surface4to2z);
title('Fit from -4 to -2 betamu');
figure(335);
surf(surface4to1y,surface4to1x,surface4to1z);
title('Fit from -4 to -1 betamu');
figure(336);
surf(surface4to0y,surface4to0x,surface4to0z);
title('Fit from -4 to -0 betamu');


%plot signatures
figure(42);
hold on;
for i=1:length(signature4to2)
    plot(signature4to2{i}.x,signature4to2{i}.sig);
end
grid on;


%kdata_880_density_SI = inputDensityArray;
%kdata_880_potential_SI = inputPotentialArray;
    
end




























%end