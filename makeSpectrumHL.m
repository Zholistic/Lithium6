function [ spectrumFunc ] = makeSpectrumHL( highIntImg, lowIntImg )
%From a high intensity image representing the 'real' atom number and a low
%intensity image both expressed in optical density (OD) create a comparison
%spectrum that is OD(H)/OD(L) vs OD(L). 

%Center the images:
toCenter = []; centeredImages = []; spectrumFunc = [];
toCenter(:,:,2) = highIntImg;
toCenter(:,:,1) = lowIntImg;
centeredImages = centerImgArray(toCenter);

highIntImg = centeredImages(:,:,2);
lowIntImg = centeredImages(:,:,1);

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




end

