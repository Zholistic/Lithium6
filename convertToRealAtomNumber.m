function [ convertedImage ] = convertToRealAtomNumber( highIntImg, lowIntImg )
%From a high intensity image representing the 'real' atom number and a low
%intensity image both expressed in optical density (OD) create a comparison
%spectrum that is OD(H)/OD(L) vs OD(L). Note the "startFit" variable which
%is the point to begin a linear fit of the spectrum. This is sensitive and
%should be begun where there is a linear spectrum.

%Evaluate the current pixels' ratio (h/l) & value of l pixel:
ratioArray = []; lowIntValueArray = [];
for i=1:length(lowIntImg(:,1))
    for j=1:length(lowIntImg(1,:))
        ratioArray(i,j) = highIntImg(i,j)/lowIntImg(i,j);
        lowIntValueArray(i,j) = lowIntImg(i,j);
    end
end

%Cap the possible ratio range from -5 to 5:
for i=1:length(ratioArray(:,1))
    for j=1:length(ratioArray(1,:))
        if ratioArray(i,j) > 5
            ratioArray(i,j) = 5;
        elseif ratioArray(i,j) < -5
            ratioArray(i,j) = -5;
        end
    end
end

%Now bin based on low intensity image pixel value:
reRatioArray = reshape(ratioArray, 1, length(ratioArray(:,1))*length(ratioArray(1,:)));
reLowIntValArray = reshape(lowIntValueArray, 1, length(lowIntValueArray(:,1))*length(lowIntValueArray(1,:)));

nbins = 100; binMean = [];
binEdges = linspace(min(min(lowIntValueArray)),max(max(lowIntValueArray)),nbins+1);

[h,whichBin] = histc(reLowIntValArray, binEdges);

for i = 1:nbins
    flagBinMembers = (whichBin == i);
    binMembers     = reRatioArray(flagBinMembers);
    binMean(i)     = mean(binMembers);
end

binLength = (max(max(lowIntValueArray))-min(min(lowIntValueArray)))/(nbins+1);

for i = 1:length(binMean)
    spectrumFunc(1,i) = binMean(i);
    spectrumFunc(2,i) = i*binLength;
end

%Now fit to the spectrum:
%Fit:
fg = @(p,x)(p(1).*x + p(2)); %function to fit with
p0 = [1 1.2];
lb = [-2 0.9];
ub = [2 2];

startFit = 40;

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
xs = spectrumFunc(2,startFit:end-1);
 coefs = lsqcurvefit(fg,p0,xs(:),spectrumFunc(1,startFit:end-1)',lb,ub,curvefitoptions);

figure(300)
plot(xs,fg(coefs,xs)); hold on; plot(spectrumFunc(2,startFit:end),spectrumFunc(1,startFit:end),'r');
hold off;

%Now multiply the pixels by the appropriate factor:

for i=1:length(lowIntImg(:,1))
    for j=1:length(lowIntImg(1,:))
       %What value to scale the pixel by?
       currPixel = lowIntImg(i,j);
       scaleFact = fg(coefs,currPixel);
       scaledImage(i,j) = currPixel*scaleFact;
    
    end
end

convertedImage = scaledImage;


end

