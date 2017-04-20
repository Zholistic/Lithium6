directory = 'C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\Datas4\';
directorylist = dir(directory);

dataNumSys = 4; %new data numbering system from set 4 on

filenumber = length(directorylist)-2;

BetaEBCells = []; TempCells = []; ChemPotCells = [];
TempArray = []; ChemPotArray = []; BetaEBArray = [];
tempCount = 1; muCount = 1; betaEBCount = 1;
for i=1:length(directorylist)-2
    
filename = directorylist(i+2).name

%BetaEBCells{i} = str2num(filename(4:5))./100;

datafilename = [directory filename];

fid = fopen(datafilename,'rt');
%C = textscan(fid, '%s', 'Delimiter','\t'); %tokenize into tab seperated tokens
if i==13 || i==12 || i==14
C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f');   
display(['hit 10 columns for' filename]);
else
C = textscan(fid,'%f %f %f %f');
end
fclose(fid);


TempCells{i} = C{1};
for k=1:length(C{1})
   TempArray(tempCount) = C{1}(k);
   tempCount = tempCount + 1;
end

ChemPotCells{i} = C{3};
for k=1:length(C{3})
   ChemPotArray(muCount) = C{3}(k);
   muCount = muCount + 1;
end

for j=1:length(TempCells{i})
    if(dataNumSys == 4)
        BetaEBCells{i}(j) = str2num(filename(6:7))./100;
        BetaEBArray(betaEBCount) = str2num(filename(6:7))./100;
        if(i > 37)
            disp('hit 100 (ie 1.0 betaeb)');
            BetaEBCells{i}(j) = str2num(filename(6:8))./100;
            BetaEBArray(betaEBCount) = str2num(filename(6:8))./100;
        end
    else
        BetaEBCells{i}(j) = str2num(filename(4:5))./100;
        BetaEBArray(betaEBCount) = str2num(filename(4:5))./100;
    end
    betaEBCount = betaEBCount +1;
end

end



%plot3(x,y,z)
close all;
for i=1:filenumber
    
hold on;
%plot3(1./(2*pi*TempCells{i}),ChemPotCells{i}./TempCells{i},BetaEBCells{i},'.');
plot3(ChemPotCells{i}./TempCells{i},BetaEBCells{i},1./(2*pi*TempCells{i}),'.');
grid on;

end
xlabel('Beta mu') % x-axis label
ylabel('Beta Eb') % y-axis label
zlabel('n (1/ (2\pi T)') % y-axis label

%-------------------------------------------
x = []; y = []; z = []; xlin = []; ylin = []; Z = []; X = []; Y = [];

x = ChemPotArray./TempArray;
y = BetaEBArray;
z = 1./(2.*pi.*TempArray);

xlin = linspace(min(x),max(x),50);
ylin = linspace(min(y),max(y),50);
        

[X,Y] = meshgrid(xlin,ylin);

Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear

surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([-11 5 0 0.5 0 2]); grid on;
xlabel('Beta mu') % x-axis label
ylabel('Beta Eb') % y-axis label
zlabel('f_n') % y-axis label
view(3);

if(0)
x = []; y = []; z = []; xlin = []; ylin = []; Z = [];
x = reshape(xdatas, 1, length(xdatas(:,1))*length(xdatas(1,:)));
y = reshape(ydatas, 1, length(ydatas(:,1))*length(ydatas(1,:)));
z = reshape(zdatas, 1, length(zdatas(:,1))*length(zdatas(1,:)));
%y = masterMatrix(2,:);
%z = masterMatrix(3,:).*2.*pixelLength./1e-6;

%reshape(logkfa2ds, 1, length(logkfa2ds(:,1))*length(logkfa2ds(1,:)));

xlin = linspace(min(x),max(x),50);
ylin = linspace(min(y),max(y),50);       

[X,Y] = meshgrid(xlin,ylin);

%Z = griddata(x,y,z,X,Y,'cubic');
Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear
end

%------------------------

betamuselect = 30;
figure(2); plot(Y(:,betamuselect),Z(:,betamuselect)); %Fn vs betaEb for a
figure(22); plot(Y(:,betamuselect),smooth(Z(:,betamuselect)));
figure(23); plot(Y(:,betamuselect),smooth(Z(:,betamuselect))-Z(:,betamuselect));
%given betamu
figure(3); plot(Y(1:49,betamuselect),diff(Z(:,betamuselect)));
figure(4); plot(Y(2:49,betamuselect),diff(Z(:,betamuselect),2));
figure(33); plot(Y(1:49,betamuselect),diff(smooth(Z(:,betamuselect))));
figure(44); plot(Y(2:49,betamuselect),diff(smooth(Z(:,betamuselect)),2));

%---------------------------------------
%Generate virial for same Beta_eb:

beta_mus = -10:0.01:5;
for i=1:filenumber
    this_betaeb = BetaEBCells{i}(1);
    b2(i) = huicoeffs(this_betaeb,1,0);
    b3(i) = huicoeffs(this_betaeb,2,0);
    
    virialtwos(i,:) = (2/(4*pi)).*(log(1+exp(beta_mus)) + 2.*b2(i).*(exp(beta_mus).^2));
    virialthrees(i,:) = (2/(4*pi)).*(log(1+exp(beta_mus)) + 2.*b2(i).*(exp(beta_mus).^2) + 3*b3(i)*(exp(beta_mus)).^3);
    betaebvirials(i,:) = ones(1,length(beta_mus)).*this_betaeb;
    
end

if(0)
betaeb = 0.2;
b2 = huicoeffs(betaeb,1,0);
b3 = huicoeffs(betaeb,2,0);
vcoefs = [b2 b3];

lambdaDB = sqrt((beta*2*pi*hbar^2)/(massL6));

%omegaz = 2*pi*5000;
%virialEqn = ['((1/(9.988e-27*' num2str(omegaz) '/ 1.05457e-34))*(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*a)))) * (log(1 + exp((((1.05457e-34*' num2str(omegaz) ')/2))*(1/((1.38e-23)*a)) * (b - x))) + (2*' num2str(vcoefs(1)) '*exp((((1.05457e-34*' num2str(omegaz) ')/2))*2*(1/((1.38e-23)*a))*(b - x))) + (3*(' num2str(vcoefs(2)) ')*exp(((1.05457e-34*' num2str(omegaz) ')/2)*3*(1/((1.38e-23)*a))*(b - x))))'];
beta_mus = -10:0.01:5;
%plot(beta_mus,(log(1+exp(beta_mus)) + 2*b2*exp(beta_mus)+3*b3*exp(beta_mus)))
figure(2); plot(ChemPotCells{10}./TempCells{10},1./(2*pi*TempCells{10}),'.');
grid on;
hold on; plot(beta_mus,(2/(4*pi)).*(log(1+exp(beta_mus)) + 2.*b2.*(exp(beta_mus).^2)));
hold on; plot(beta_mus,(2/(4*pi)).*(log(1+exp(beta_mus)) + 2.*b2.*(exp(beta_mus).^2) + 3*b3*(exp(beta_mus)).^3));
axis([-5 0 0 0.25]);
title(['F_n vs BetaMu for the Viral, GG at 0.2 BetaEb']);

n_virial = (1/(lambdaDB^2))*(log(1+exp(beta_mus)) + 2*b2*(exp(beta_mus)).^2 + 3*b3*(exp(beta_mus)).^3);
end

%plot3(x,y,z)
close all;
for i=1:filenumber
    
hold on;
%plot3(1./(2*pi*TempCells{i}),ChemPotCells{i}./TempCells{i},BetaEBCells{i},'.');
plot3(ChemPotCells{i}./TempCells{i},BetaEBCells{i},1./(2*pi*TempCells{i}),'.');
plot3(beta_mus,betaebvirials(i,:),virialtwos(i,:));
plot3(beta_mus,betaebvirials(i,:),virialthrees(i,:));
grid on;

end
axis([-11 5 0 0.5 0 2]); grid on;
xlabel('Beta mu') % x-axis label
ylabel('Beta Eb') % y-axis label
zlabel('F_n') % y-axis label
view(3);

betaMusGG = []; fnGG = [];
for i=1:filenumber
betaMusGG{i} = ChemPotCells{i}./TempCells{i};
fnGG{i} = 1./(2*pi*TempCells{i});
end

composite_TfuncFn = []; toChooseBetaMusVirial = []; toChooseBetaMusGG = [];
composite_TfuncGG = []; composite_TfuncBetaMu = [];
for i=1:filenumber
    
    %lhs of function:
    %first find element of the virial closest to -3.5 BetaMu
    toChooseBetaMusVirial{i} = find(beta_mus <= -3.5);
    composite_TfuncVirial{i} = virialthrees(i,toChooseBetaMusVirial{i});
    
    %rhs of function:
    %find elements of the GG > -3.5 BetaMu
    %toChooseBetaMusGG{i} = find(betaMusGG{i} > -3.5);
    toChooseBetaMusGG{i} = find(betaMusGG{i} > -3.5 & betaMusGG{i} < 1);
    composite_TfuncGG{i} = fnGG{i}(toChooseBetaMusGG{i});
    
    composite_TfuncFn{i} = cat(1,composite_TfuncVirial{i}',flip(composite_TfuncGG{i}));
    composite_TfuncBetaMu{i} = cat(1,beta_mus(toChooseBetaMusVirial{i})',flip(betaMusGG{i}(toChooseBetaMusGG{i})));
    
    if(~issorted(composite_TfuncFn{i}))
        composite_TfuncFn{i} = sort(composite_TfuncFn{i});
        composite_TfuncBetaMu{i} = sort(composite_TfuncBetaMu{i});
    end
       
end

figure(3000)
for i=1:filenumber
    hold on;
    plot(composite_TfuncBetaMu{i},composite_TfuncFn{i},'.');
    %Fnsorted = issorted(composite_TfuncFn{i}')
    %BetaMusorted = issorted(composite_TfuncBetaMu{i})
end
grid on;
    


%-----------
%Fit function to composite points:
%Function of form (a + a2x + a3x^2 ...)/(1 + e^(-bx + c))
%(a + b*x + c*x^2 + d*x^3 + e*x^4 + f*x^5 + g*x^6)/(1 + exp(-l*x + m))

for i=1:filenumber

x = composite_TfuncBetaMu{i};
y = composite_TfuncFn{i};

p0 = [0.365 0.0986 0.013 -0.0002 -7e-05 5.9e-06 -1.27e-07 0.8922 0.4494];
lb = [0.025 0.001 0.001 -0.01 -1e-03 1e-08 -1e-03 0.01 0.01];
ub = [1 1 0.1 0.001 1e03 1e03 1e02 2 2];

fitEqn = ['(a + b*x + c*x^2 + d*x^3 + e*x^4 + f*x^5 + g*x^6)/(1 + exp(-l*x + m))'];
[f1, gof, foutput] = fit(x,y,fitEqn,'Start', p0,'Lower',lb,'Upper',ub,'MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-11,'TolX',1e-11);

fitFuncsComp{i} = f1;
gofFuncsComp{i} = gof;
foutputFuncsComp{i} = foutput;



end

close all;
figure(4000)
axis([-11 15 0 3]); grid on;
for i=1:filenumber
    hold on;
    plot(fitFuncsComp{i},'.');
    %Fnsorted = issorted(composite_TfuncFn{i}')
    %BetaMusorted = issorted(composite_TfuncBetaMu{i})
end

%plot3(x,y,z)
compxs = -10:0.1:1;
figure(4001);
close all;
for i=1:filenumber
    
hold on;
%plot3(1./(2*pi*TempCells{i}),ChemPotCells{i}./TempCells{i},BetaEBCells{i},'.');
plot3(ChemPotCells{i}./TempCells{i},BetaEBCells{i},1./(2*pi*TempCells{i}),'.');
plot3(beta_mus,betaebvirials(i,:),virialtwos(i,:));
plot3(beta_mus,betaebvirials(i,:),virialthrees(i,:));
plot3(compxs,ones(length(compxs),1).*betaebvirials(i,1),feval(fitFuncsComp{i}, compxs));
grid on;

end
axis([-11 15 0 1 0 2]); grid on;
xlabel('Beta mu') % x-axis label
ylabel('Beta Eb') % y-axis label
zlabel('F_n') % y-axis label
view(3);

    
    
    
    
    
    
    

